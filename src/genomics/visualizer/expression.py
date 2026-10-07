"""How much of each gene product an individual makes: relative to the reference genome and absolute.

Two different questions with different evidence:

**Relative** (fold change against the reference genome). AlphaGenome predicts RNA-seq coverage of
the reference window and of each of the individual's haplotypes with the same model and tissue
track; their ratio over the gene's exons is the predicted change in mRNA caused by the individual's
sequence. Ratios cancel most of the model's calibration error (e.g. the CAGE/RNA bias between
tissues), so this is the quantity the model can be trusted with. The individual (diploid) fold is
``(H1 + H2) / (2 x reference)``; per-haplotype folds show allele-specific expression.

**Absolute** (TPM-like units). AlphaGenome's coverage is not a concentration, so an absolute level
needs an observed anchor: the tissue's measured level of the gene in a reference population (GTEx v8
median TPM for lymphoblastoid cells, Human Protein Atlas single-cell nCPM for melanocytes). The
reference genome is assumed to sit at that level and the individual at ``anchor x fold``. It is an
estimate for this genome, not a measurement of this person.

**Protein.** No model here predicts protein abundance from sequence, so protein is *relative only*:
the mRNA fold times the change in the share of transcripts expected to make a protein from the
annotated start without nonsense-mediated decay (NMD-targeted and start-lost products, including
candidate isoforms weighted by their predicted junction usage, make none). It assumes unchanged
translation efficiency and protein stability; a missense change alters which protein is made, not
how much.

Predictions come from the dataset's stored ``rna_seq`` when the reference and both haplotypes have
tracks for a tissue, otherwise from AlphaGenome on demand (cached under ``<cache>/expression``).
"""
from __future__ import annotations

import csv
import io
import json
import re
import threading
import zipfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from genomics.visualizer.alphagenome import AlphaGenomeUnavailable, predict_cached
from genomics.visualizer.datasets import Dataset
from genomics.visualizer.products import HAPLOTYPES, REFERENCE
from genomics.visualizer.remote import DAY, RemoteCache
from genomics.workflows.alphagenome.outputs import SUPPORTED_LENGTHS

ProgressFn = Callable[[float, str], None]

HPA_SINGLE_CELL_URL = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip"
HPA_CELL_LINE_URL = "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip"
# kind -> (download, column holding the cell type / line, value column, unit, label)
HPA_TABLES = {
    "single_cell": (HPA_SINGLE_CELL_URL, 2, 3, "nCPM", "Human Protein Atlas single-cell RNA, {name} (median of donors)"),
    "cell_line": (HPA_CELL_LINE_URL, 2, 5, "nTPM", "Human Protein Atlas cell line RNA, {name}"),
}
GTEX_LCL = "Cells_EBV-transformed_lymphocytes"
GTEX_SKIN = "Skin_Sun_Exposed_Lower_leg"
HPA_ANCHOR = {"source": "hpa_single_cell", "cell_type": "melanocytes", "unit": "nCPM",
              "label": "Human Protein Atlas single-cell RNA, melanocytes (median of donors)"}
GTEX_ANCHOR = {"source": "gtex", "tissue": GTEX_LCL, "unit": "TPM", "label": "GTEx v8 median, EBV-transformed lymphocytes"}
SKIN_SHARES = {"tissue": GTEX_SKIN, "label": "GTEx v8 median transcript TPM, sun-exposed skin (bulk tissue standing in for melanocytes)", "proxy": True}
LCL_SHARES = {"tissue": GTEX_LCL, "label": "GTEx v8 median transcript TPM, EBV-transformed lymphocytes", "proxy": False}
POOLED_SHARES = {"tissue": None, "label": "GTEx v8 median transcript TPM averaged over all GTEx tissues (no tissue-matched transcript data)", "proxy": True}

# The tissues shown by default (curated anchors); any other AlphaGenome ontology term is resolved by
# :func:`tissue_spec` from the track catalog, GTEx's tissue list and the Human Protein Atlas.
TISSUES: List[Dict[str, Any]] = [
    {"ontology": "CL:1000458", "label": "melanocyte of skin", "group": "melanocyte", "anchor": HPA_ANCHOR, "transcripts": SKIN_SHARES},
    {"ontology": "CL:2000045", "label": "foreskin melanocyte", "group": "melanocyte", "anchor": HPA_ANCHOR, "transcripts": SKIN_SHARES},
    {"ontology": "EFO:0000572", "label": "lymphoblastoid cells (GTEx)", "group": "LCL", "anchor": GTEX_ANCHOR, "transcripts": LCL_SHARES},
    {"ontology": "EFO:0002784", "label": "GM12878 (LCL)", "group": "LCL", "anchor": {**GTEX_ANCHOR, "proxy": True}, "transcripts": LCL_SHARES},
]
DEFAULT_ONTOLOGIES = [t["ontology"] for t in TISSUES]
SIMILAR = 0.10  # within ±10% of the reference counts as "about the same"
EXPRESSED_ANCHOR = 1.0  # TPM / nCPM below which the tissue is taken not to express the gene
MAX_NMD_SHARE = 0.95
PROTEIN_ABSOLUTE_NOTE = ("No absolute protein level: there is no tissue-matched quantitative proteomics reference wired in, and "
                         "protein per mRNA varies by orders of magnitude between genes. Only the relative change is shown.")


def cell_type_key(name: str) -> str:
    """'melanocyte of skin' / 'Melanocytes' -> 'melanocyte'; 'T-cell' / 't-cells' -> 't-cell'."""
    s = re.sub(r"\s+of\s+.*$", "", str(name).strip().lower())
    s = re.sub(r"\s+", " ", s)
    return s[:-1] if s.endswith("s") and not s.endswith("ss") else s


def cell_line_key(name: str) -> str:
    """'Hep-G2' / 'HepG2' -> 'hepg2'; 'K-562' / 'K562' -> 'k562'."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


class ReferenceLevels:
    """Observed expression of a gene in a tissue, used as the reference genome's absolute level."""

    def __init__(self, remote: RemoteCache, gtex, cache_dir: Optional[Path]):
        self.remote = remote
        self.gtex = gtex
        self.cache_dir = Path(cache_dir) / "expression" if cache_dir else None
        self._columns: Dict[Tuple[str, str], Dict[str, float]] = {}
        self._names: Dict[str, List[str]] = {}
        self._lock = threading.Lock()

    def _path(self, name: str) -> Optional[Path]:
        return self.cache_dir / name if self.cache_dir else None

    def _scan(self, kind: str, wanted: Optional[str]) -> None:
        """One pass over an HPA table: its column names, and the values of column ``wanted``."""
        url, name_col, value_col, _unit, _label = HPA_TABLES[kind]
        raw = self.remote.get_bytes(url, ttl=180 * DAY)
        names = set()
        values: Dict[str, float] = {}
        with zipfile.ZipFile(io.BytesIO(raw)) as zf:
            with zf.open(zf.namelist()[0]) as fh:
                reader = csv.reader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t")
                next(reader, None)
                for row in reader:
                    if len(row) <= value_col:
                        continue
                    names.add(row[name_col])
                    if row[name_col] == wanted:
                        try:
                            values[row[0]] = float(row[value_col])
                        except ValueError:
                            pass
        self._names[kind] = sorted(names)
        path = self._path(f"hpa_{kind}_names.json")
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(self._names[kind]), encoding="utf-8")
        if wanted is not None:
            self._store(kind, wanted, values)

    def _store(self, kind: str, column: str, values: Dict[str, float]) -> None:
        self._columns[(kind, column)] = values
        path = self._path(f"hpa_{kind}_{re.sub(r'[^A-Za-z0-9]+', '_', column)}.json")
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(values), encoding="utf-8")

    def names(self, kind: str) -> List[str]:
        """Cell types (``single_cell``) or cell lines (``cell_line``) of an HPA table."""
        with self._lock:
            if kind not in self._names:
                path = self._path(f"hpa_{kind}_names.json")
                if path is not None and path.exists():
                    self._names[kind] = json.loads(path.read_text(encoding="utf-8"))
                else:
                    self._scan(kind, None)
            return self._names[kind]

    def column(self, kind: str, column: str) -> Dict[str, float]:
        """{Ensembl gene id: value} of one HPA cell type / cell line (one pass over the table, then cached)."""
        with self._lock:
            hit = self._columns.get((kind, column))
            if hit is not None:
                return hit
            path = self._path(f"hpa_{kind}_{re.sub(r'[^A-Za-z0-9]+', '_', column)}.json")
            if path is not None and path.exists():
                self._columns[(kind, column)] = json.loads(path.read_text(encoding="utf-8"))
                return self._columns[(kind, column)]
            self._scan(kind, column)
            return self._columns.get((kind, column), {})

    def _hpa_index(self, cell_type: str) -> Dict[str, float]:
        return self.column("single_cell", cell_type)

    def match_hpa(self, name: str, biosample_type: str) -> Optional[Dict[str, Any]]:
        """An HPA anchor for an AlphaGenome biosample, matched by name (cell lines, else cell types)."""
        kinds = ("cell_line",) if biosample_type == "cell_line" else ("single_cell",) if biosample_type in ("primary_cell", "in_vitro_differentiated_cells") else ()
        for kind in kinds:
            key = cell_line_key if kind == "cell_line" else cell_type_key
            target = key(name)
            match = next((n for n in self.names(kind) if key(n) == target), None)
            if match:
                _url, _nc, _vc, unit, label = HPA_TABLES[kind]
                field = "cell_line" if kind == "cell_line" else "cell_type"
                return {"source": f"hpa_{kind}", field: match, "unit": unit, "label": label.format(name=match), "matched_by_name": True}
        return None

    def level(self, anchor: Optional[Dict[str, Any]], symbol: str, gene_id: Optional[str]) -> Dict[str, Any]:
        if not anchor:
            return {"unit": None, "source": None, "proxy": False, "value": None,
                    "error": "no observed level for this tissue (not a GTEx tissue, and no Human Protein Atlas cell type or cell line of that name)"}
        out = {"unit": anchor["unit"], "source": anchor["label"], "proxy": bool(anchor.get("proxy")), "value": None,
               "matched_by_name": bool(anchor.get("matched_by_name"))}
        try:
            if anchor["source"] in ("hpa_single_cell", "hpa_cell_line"):
                if not gene_id:
                    out["error"] = "no Ensembl gene id in the annotation"
                    return out
                kind = anchor["source"][4:]
                value = self.column(kind, anchor.get("cell_type") or anchor.get("cell_line")).get(gene_id.split(".")[0])
                out["value"] = value
                if value is None:
                    out["error"] = "gene not in the Human Protein Atlas table"
            elif anchor["source"] == "gtex":
                gene = self.gtex.gene(symbol)
                # only an exact symbol: a near match would anchor this gene to another gene's level
                if not gene or not gene.get("gencodeId") or str(gene.get("geneSymbol", "")).upper() != symbol.upper():
                    out["error"] = f"{symbol} not found in GTEx"
                    return out
                out["value"] = self.gtex.median_expression(gene["gencodeId"], anchor["tissue"])
        except Exception as exc:  # a missing anchor (offline, service down) never fails the page
            out["error"] = str(exc)
        return out


def tissue_spec(curie: str, info: Optional[Dict[str, Any]], references: Optional["ReferenceLevels"] = None,
                gtex_by_ontology: Optional[Dict[str, str]] = None, gtex_by_name: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """How the products pipeline treats an AlphaGenome ontology term.

    ``info`` (from the track catalog): ``name``, ``type``, ``gtex_tissue``. The absolute anchor and
    the transcript shares come from GTEx when the term is a GTEx tissue (its RNA-seq tracks are GTEx's,
    GTEx lists the same ontology term, or a tissue has the name of a single-site GTEx tissue, e.g.
    ENCODE's "liver" and GTEx's "Liver"); otherwise the anchor is a Human Protein Atlas cell type or
    cell line of the same name (resolved with ``references``, which may download the table), and the
    transcript shares are GTEx's average over all tissues.
    """
    known = next((t for t in TISSUES if t["ontology"] == curie), None)
    if known is not None:
        return known
    info = info or {}
    name = info.get("name") or curie
    gtex = info.get("gtex_tissue") or (gtex_by_ontology or {}).get(curie)
    by_name = False
    if not gtex and info.get("type") == "tissue" and gtex_by_name:
        from genomics.visualizer.gtex import tissue_name_key

        gtex = gtex_by_name.get(tissue_name_key(name))
        by_name = bool(gtex)
    spec: Dict[str, Any] = {"ontology": curie, "label": name, "group": (info.get("type") or "other").replace("_", " "), "type": info.get("type")}
    if gtex:
        spec["anchor"] = {"source": "gtex", "tissue": gtex, "unit": "TPM", "label": f"GTEx v8 median, {gtex.replace('_', ' ')}", "matched_by_name": by_name}
        spec["transcripts"] = {"tissue": gtex, "label": f"GTEx v8 median transcript TPM, {gtex.replace('_', ' ')}", "proxy": False}
        spec["gtex_tissue"] = gtex
    else:
        anchor = None
        if references is not None:
            try:
                anchor = references.match_hpa(name, info.get("type") or "")
            except Exception:  # HPA unreachable: relative only
                anchor = None
        spec["anchor"] = anchor
        spec["transcripts"] = POOLED_SHARES
    return spec


def _tracks_for(records: List[Dict[str, Any]], ontology: str, strand: str) -> List[int]:
    return [i for i, r in enumerate(records) if r.get("ontology_curie") == ontology and str(r.get("strand", ".")) in (strand, ".")]


def _read_meta(npz: Path) -> List[Dict[str, Any]]:
    meta = npz.with_name(f"{npz.stem}_metadata.json")
    if not meta.exists():
        return []
    with open(meta, "r", encoding="utf-8") as f:
        payload = json.load(f)
    return list(payload.get("metadata") or []) if isinstance(payload, dict) else list(payload)


def exon_offsets(transcripts: List[Dict[str, Any]], length: int) -> np.ndarray:
    """Reference offsets covered by any exon of the gene's (complete) transcripts."""
    mask = np.zeros(length, dtype=bool)
    for t in transcripts:
        for s, e in t["exons"]:
            mask[max(int(s), 0):min(int(e), length)] = True
    return np.nonzero(mask)[0]


def mean_over(values: np.ndarray, local: np.ndarray, offsets: np.ndarray, columns: List[int]) -> float:
    """Mean predicted coverage over ``offsets`` (reference) read at the haplotype's own positions."""
    idx = local[offsets]
    idx = idx[(idx >= 0) & (idx < values.shape[0])]
    if not idx.size or not columns:
        return float("nan")
    return float(values[np.ix_(idx, columns)].mean())


def fold_label(fold: Optional[float]) -> str:
    if fold is None or not np.isfinite(fold):
        return "unknown"
    if fold >= 1 + SIMILAR:
        return "higher"
    if fold <= 1 / (1 + SIMILAR):
        return "lower"
    return "similar"


def _ratio(a: float, b: float) -> Optional[float]:
    return float(a / b) if np.isfinite(a) and np.isfinite(b) and b > 0 else None


def productive_share(product: Dict[str, Any], candidates: List[Dict[str, Any]], hap: str, column: Optional[int]) -> float:
    """Share of the gene's transcripts expected to make a protein (annotated start, no NMD)."""
    if not product.get("coding") or not product.get("protein") or (product.get("nmd") or {}).get("predicted"):
        return 0.0
    lost = 0.0
    if column is not None:
        for c in candidates:
            p = (c.get("products") or {}).get(hap) or {}
            junction = (c.get("junction") or {}).get(hap) or (c.get("junction") or {}).get(REFERENCE) or {}
            usage = (junction.get("usage") or [0.0] * (column + 1))[column]
            if (p.get("nmd") or {}).get("predicted") or not p.get("protein"):
                lost += float(usage)
    return float(1.0 - min(lost, MAX_NMD_SHARE))


class ExpressionService:
    def __init__(self, app):
        self.app = app
        self.references = ReferenceLevels(app.remote, app.gtex, app.cache_dir)
        self._memory: Dict[Tuple, Dict[str, Any]] = {}
        self._junctions: Dict[Tuple, Dict[str, Any]] = {}

    @staticmethod
    def _ontologies(tissues: Optional[List[Dict[str, Any]]]) -> Tuple[str, ...]:
        return tuple(t["ontology"] for t in (tissues or TISSUES))

    def key(self, dataset: Dataset, sample: str, gene: str, target: Optional[str], tissues: Optional[List[Dict[str, Any]]] = None) -> str:
        return f"expression:{dataset.id}:{sample}:{gene}:{target or ''}:{','.join(self._ontologies(tissues))}"

    def cached(self, dataset: Dataset, sample: str, gene: str, target: Optional[str], tissues: Optional[List[Dict[str, Any]]] = None) -> Optional[Dict[str, Any]]:
        return self._memory.get((dataset.fingerprint, sample, gene, target, self._ontologies(tissues)))

    def junctions(self, dataset: Dataset, sample: str, gene: str, ontology: str) -> Optional[Dict[str, Any]]:
        """Splice junctions predicted on demand with a tissue's RNA-seq: {hap: (arrays, track records)}."""
        return self._junctions.get((dataset.fingerprint, sample, gene, ontology))

    def _predictions(self, dataset: Dataset, sample: str, gene: str, strand: str, views, progress: ProgressFn,
                     specs: Optional[List[Dict[str, Any]]] = None) -> Tuple[Dict[str, Dict[str, Any]], List[str]]:
        """{ontology: {"source", "values": {hap: (matrix, columns)}}} and per-tissue problems.

        Tissues missing from the stored predictions are predicted one at a time (RNA-seq and splice
        junctions in the same call, so a tissue's cache entry is reused whatever else is asked for)."""
        stored_paths = {REFERENCE: dataset.reference_prediction_path(gene, "rna_seq")}
        stored_paths.update({hap: dataset.prediction_path(sample, gene, hap, "rna_seq") for hap in HAPLOTYPES})
        stored_meta = {h: _read_meta(p) if p.exists() else [] for h, p in stored_paths.items()}
        out: Dict[str, Dict[str, Any]] = {}
        missing: List[str] = []
        for tissue in specs or TISSUES:
            cols = {h: _tracks_for(stored_meta[h], tissue["ontology"], strand) for h in stored_paths}
            if all(cols.values()):
                out[tissue["ontology"]] = {"source": "stored", "columns": cols}
            else:
                missing.append(tissue["ontology"])
        problems: List[str] = []
        loaded: Dict[str, np.ndarray] = {}
        for ontology, entry in out.items():
            entry["values"] = {}
            for h, path in stored_paths.items():
                if h not in loaded:
                    from genomics.visualizer.signals import load_prediction_matrix

                    loaded[h] = load_prediction_matrix(path)
                entry["values"][h] = (loaded[h], entry["columns"][h])
        if missing:
            cache_dir = Path(self.app.cache_dir) / "expression" if self.app.cache_dir else None
            sequences = {REFERENCE: self.app.signals.reference_sequence(dataset, gene)}
            sequences.update({h: self.app.signals.haplotype_sequence(dataset, sample, gene, h) for h in HAPLOTYPES})
            bad = sorted({len(seq) for seq in sequences.values()} - set(SUPPORTED_LENGTHS))
            if bad:
                problems.append(f"{', '.join(missing)}: not in the stored predictions, and a {bad[0]:,}-bp window is not an AlphaGenome input length")
                missing = []
            for n, ontology in enumerate(missing):
                predicted: Dict[str, Any] = {}
                try:
                    for i, (h, seq) in enumerate(sequences.items()):
                        progress(0.05 + 0.75 * (n * 3 + i) / (3 * len(missing)),
                                 f"AlphaGenome RNA-seq and splice junctions of the {'reference' if h == REFERENCE else h} window in {ontology}")
                        predicted[h] = predict_cached(getattr(self.app, "alphagenome", None), cache_dir, seq, ["rna_seq", "splice_junctions"], [ontology])
                except AlphaGenomeUnavailable as exc:
                    problems.append(f"{ontology}: not in the stored predictions and AlphaGenome is unavailable ({exc})")
                    break
                except Exception as exc:  # a failed prediction is reported per tissue, not as a page error
                    problems.append(f"{ontology}: AlphaGenome prediction failed ({type(exc).__name__}: {exc})")
                    continue
                cols = {h: _tracks_for(p["rna_seq"][1], ontology, strand) for h, p in predicted.items()}
                if all(cols.values()):
                    out[ontology] = {"source": "on_demand", "columns": cols, "values": {h: (predicted[h]["rna_seq"][0], cols[h]) for h in predicted}}
                else:
                    problems.append(f"{ontology}: AlphaGenome has no RNA-seq track for this tissue on the {strand} strand")
                self._junctions[(dataset.fingerprint, sample, gene, ontology)] = {h: p["splice_junctions"] for h, p in predicted.items()}
        return out, problems

    def compute(self, dataset: Dataset, sample: str, gene: str, models: List[Dict[str, Any]], target: Optional[str], progress: ProgressFn,
                tissues: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
        specs = tissues or TISSUES
        products = self.app.products
        progress(0.02, "Gene products")
        payload = products.products(dataset, sample, gene, models, target=target)
        if not payload.get("transcripts"):
            return {"available": False, "message": payload.get("message") or "No transcripts"}
        chosen = next((g for g in models if g["name"] == payload["target"]), {})
        views = products.views(dataset, sample, gene)
        strand = payload["strand"]
        length = views[REFERENCE].local.size
        offsets = exon_offsets(payload["transcripts"], length)
        predictions, problems = self._predictions(dataset, sample, gene, strand, views, progress, specs)
        splicing = payload.get("splicing") or {}
        junction_tracks = splicing.get("tracks") or []
        base = next(t for t in payload["transcripts"] if t["id"] == payload["base"])
        candidates = splicing.get("candidates") or []
        tissues = []
        for i, tissue in enumerate(specs):
            progress(0.85 + 0.1 * i / max(len(specs), 1), f"Reference level in {tissue['label']}")
            row: Dict[str, Any] = {"ontology": tissue["ontology"], "label": tissue["label"], "group": tissue["group"]}
            entry = predictions.get(tissue["ontology"])
            anchor = self.references.level(tissue.get("anchor"), payload["target"], chosen.get("id"))
            row["anchor"] = anchor
            if entry is None:
                row["available"] = False
                tissues.append(row)
                continue
            row["available"] = True
            row["source"] = entry["source"]
            coverage = {}
            for h in (REFERENCE,) + HAPLOTYPES:
                values, cols = entry["values"][h]
                coverage[h] = mean_over(values, views[h].local, offsets, cols)
            fold = {h: _ratio(coverage[h], coverage[REFERENCE]) for h in HAPLOTYPES}
            individual = _ratio(coverage["H1"] + coverage["H2"], 2 * coverage[REFERENCE])
            row["mrna"] = {
                "predicted": {h: round(v, 4) if np.isfinite(v) else None for h, v in coverage.items()},
                "fold": {**fold, "individual": individual},
                "verdict": fold_label(individual),
            }
            a = anchor.get("value")
            row["expressed"] = None if a is None else bool(a >= EXPRESSED_ANCHOR)
            if a is not None:
                row["mrna"]["absolute"] = {
                    REFERENCE: a,
                    **{h: (a / 2 * fold[h]) if fold[h] is not None else None for h in HAPLOTYPES},
                    "reference_haplotype": a / 2,
                    "individual": a * individual if individual is not None else None,
                }
            # isoform shares from the junction tracks of the same tissue
            column = next((t["index"] for t in junction_tracks if t.get("ontology") == tissue["ontology"]), None)
            row["junction_track"] = column
            prod = {h: productive_share((base["products"].get(h) or {}), candidates, h, column) for h in (REFERENCE,) + HAPLOTYPES}
            protein = {h: (fold[h] * prod[h] / prod[REFERENCE]) if fold[h] is not None and prod[REFERENCE] > 0 else None for h in HAPLOTYPES}
            row["protein"] = {
                "productive_share": {h: round(v, 3) for h, v in prod.items()},
                "fold": {**protein, "individual": (protein["H1"] + protein["H2"]) / 2 if None not in protein.values() else None},
                "variant": {h: (base["products"].get(h) or {}).get("change") for h in HAPLOTYPES},
                "absolute": None,
            }
            row["protein"]["verdict"] = fold_label(row["protein"]["fold"]["individual"])
            row["isoforms"] = _isoform_rows(base, candidates, column)
            tissues.append(row)
        result = {
            "available": True, "sample": sample, "gene": gene, "target": payload["target"], "base_name": payload["base_name"],
            "tissues": tissues, "problems": problems, "exonic_bases": int(offsets.size),
            "notes": {
                "relative": "Fold change against the reference genome, predicted by AlphaGenome for the same tissue track (mean RNA-seq coverage over the gene's exons).",
                "absolute": "Relative fold × the tissue's observed level of the gene in a reference population (the reference genome is assumed to sit at that level). An estimate for this genome, not a measurement of this person.",
                "protein": "Protein (relative) = mRNA fold × change in the share of transcripts that make a protein (NMD-targeted and start-lost products make none). Assumes unchanged translation and protein stability.",
                "protein_absolute": PROTEIN_ABSOLUTE_NOTE,
                "similar": SIMILAR,
            },
        }
        self._memory[(dataset.fingerprint, sample, gene, target, self._ontologies(specs))] = result
        progress(1.0, "Expression ready")
        return result


def _isoform_rows(base: Dict[str, Any], candidates: List[Dict[str, Any]], column: Optional[int]) -> List[Dict[str, Any]]:
    """Base transcript and candidate isoforms with their junction usage on this tissue (shares)."""
    if column is None:
        return []
    rows = []
    for c in candidates:
        j = c.get("junction") or {}
        usage = {h: ((j.get(h) or {}).get("usage") or [None] * (column + 1))[column] for h in (REFERENCE,) + HAPLOTYPES}
        if max(v or 0 for v in usage.values()) < 0.05:
            continue
        rows.append({"kind": c.get("kind"), "start": j.get("start"), "end": j.get("end"), "share": usage,
                     "nmd": {h: bool(((c.get("products") or {}).get(h) or {}).get("nmd", {}).get("predicted")) for h in (REFERENCE,) + HAPLOTYPES}})
    rows.sort(key=lambda r: -max(v or 0 for v in r["share"].values()))
    return rows[:4]

