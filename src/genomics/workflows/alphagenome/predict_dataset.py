"""AlphaGenome predictions for every haplotype window of a canonical-layout dataset.

    python -m genomics.workflows.alphagenome.predict_dataset DATASET_DIR \
        --outputs RNA_SEQ,CAGE --ontology CL:1000458,UBERON:0002107 [--genes OCA2,HERC2] \
        [--samples-file ids.txt] [--haplotypes H1,H2,ref] [--overwrite]

Every AlphaGenome output can be requested (``--outputs all``): per-base tracks, 128 bp ChIP-seq
tracks, splice junctions and contact maps (see :mod:`genomics.workflows.alphagenome.outputs` for
how each is stored). Reads ``individuals/<sample>/windows/<gene>/<sample>.<H>.window.fixed.fa`` and
writes ``predictions_<H>/<output>.npz`` plus ``<output>_metadata.json`` (``{"metadata": [one record
per track], ...}``), the files the 1000 Genomes builder produces. The haplotype ``ref`` predicts each
window's reference sequence (``references/windows/<gene>/ref.window.fa`` ->
``references/windows/<gene>/predictions_ref/``, same files), the baseline the visualizer draws next
to individuals and observed data. Windows that already have every requested output are skipped,
so an interrupted run resumes where it stopped. The backend is the hosted API (``ALPHAGENOME_API_KEY``) or the self-hosted server named by
``ALPHAGENOME_ADDRESS`` (see :mod:`genomics.core.alphagenome_connection`).

Progress is printed as ``@@progress <fraction> <message>`` lines for the visualizer's task runner.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.workflows.alphagenome.outputs import ALL_OUTPUTS, OUTPUT_SPECS, SUPPORTED_LENGTHS, normalize_outputs

REFERENCE = "ref"
DETAIL_KEYS = ("biosample_name", "biosample_type", "biosample_life_stage", "gtex_tissue", "data_source", "endedness", "genetically_modified")


def emit_progress(fraction: float, message: str) -> None:
    print(f"@@progress {max(0.0, min(1.0, fraction)):.4f} {message}", flush=True)


def read_fasta(path: Path) -> str:
    with open(path, "r", encoding="utf-8") as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">")).upper()


def _clean(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def metadata_records(track_data: Any) -> List[Dict[str, Any]]:
    metadata = getattr(track_data, "metadata", None)
    if metadata is None:
        return []
    try:
        records = metadata.to_dict(orient="records")
    except AttributeError:
        return []
    return [{str(k): _clean(v) for k, v in record.items()} for record in records]


def junction_arrays(junctions: Any) -> Dict[str, np.ndarray]:
    """``JunctionData`` as arrays: 0-based half-open window offsets, strands and (junctions, tracks) values."""
    raw = getattr(junctions, "junctions", None)  # numpy array of Interval objects
    items = [] if raw is None else list(raw)
    return {
        "starts": np.asarray([int(j.start) for j in items], dtype=np.int64),
        "ends": np.asarray([int(j.end) for j in items], dtype=np.int64),
        "strands": np.asarray([str(getattr(j, "strand", ".")) for j in items], dtype="<U1"),
        "values": np.asarray(junctions.values, dtype=np.float32).reshape(len(items), -1),
    }


def write_atomic_json(path: Path, payload: Any) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


class DatasetPredictor:
    def __init__(
        self,
        dataset_dir: Path,
        outputs: Sequence[str],
        ontology_terms: Optional[Sequence[str]],
        genes: Optional[Sequence[str]] = None,
        samples: Optional[Sequence[str]] = None,
        haplotypes: Sequence[str] = ("H1", "H2"),
        overwrite: bool = False,
        timeout: float = 600.0,
        max_attempts: int = 3,
        rate_limit_delay: float = 0.0,
    ):
        self.dataset_dir = Path(dataset_dir).expanduser().resolve()
        self.meta_path = self.dataset_dir / "dataset_metadata.json"
        if not self.meta_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json not found in {self.dataset_dir}")
        self.metadata = json.loads(self.meta_path.read_text(encoding="utf-8"))
        self.outputs = normalize_outputs(outputs)
        if not self.outputs:
            raise ValueError(f"Choose outputs among {', '.join(ALL_OUTPUTS)} (or 'all')")
        self.ontology_terms = [t.strip() for t in (ontology_terms or []) if t.strip()] or None
        # Windows added after the dataset was built (e.g. by ``vcf-import --extend``) may be missing
        # from the metadata's gene list, so the window directories count too.
        windows_dir = self.dataset_dir / "references" / "windows"
        on_disk = sorted(p.name for p in windows_dir.iterdir() if p.is_dir()) if windows_dir.is_dir() else []
        listed = [str(g) for g in self.metadata.get("genes") or []]
        all_genes = listed + sorted(set(on_disk) - set(listed))
        self.genes = [g for g in (genes or all_genes) if g in all_genes]
        missing_genes = [g for g in (genes or []) if g not in all_genes]
        if missing_genes:
            raise ValueError(f"Unknown windows: {', '.join(missing_genes)}")
        all_samples = [str(s if not isinstance(s, dict) else s.get("sample_id")) for s in self.metadata.get("individuals") or []]
        self.samples = [s for s in (samples or all_samples) if s in set(all_samples)]
        self.haplotypes = [h for h in haplotypes if h in ("H1", "H2")]
        self.reference = any(str(h).lower() in (REFERENCE, "reference") for h in haplotypes)
        self.overwrite = overwrite
        self.timeout = timeout
        self.max_attempts = max_attempts
        self.rate_limit_delay = rate_limit_delay
        self.details: Dict[str, Dict[str, Any]] = {}
        self.done_outputs: Dict[str, set] = {}
        self.reference_outputs: Dict[str, set] = {}

    def _pending(self, sample: Optional[str], gene: str, hap: str) -> Tuple[Optional[Path], List[str]]:
        if hap == REFERENCE:
            case = self.dataset_dir / "references" / "windows" / gene
            fasta = case / "ref.window.fa"
        else:
            case = self.dataset_dir / "individuals" / str(sample) / "windows" / gene
            fasta = case / f"{sample}.{hap}.window.fixed.fa"
        if not fasta.exists():
            return None, []
        pred_dir = case / f"predictions_{hap}"
        if self.overwrite:
            return fasta, list(self.outputs)
        return fasta, [o for o in self.outputs if not (pred_dir / f"{o.lower()}.npz").exists()]

    def _client(self):
        from genomics.core.alphagenome_connection import create_dna_client

        return create_dna_client(timeout=60.0)

    def _save(self, case: Path, hap: str, prediction: Any, outputs: List[str], gene: str) -> List[str]:
        pred_dir = case / f"predictions_{hap}"
        pred_dir.mkdir(exist_ok=True)
        saved = []
        for name in outputs:
            spec = OUTPUT_SPECS[name]
            data = getattr(prediction, spec.attr, None)
            if data is None or getattr(data, "values", None) is None:
                print(f"[WARN] AlphaGenome returned no {name} for {case.name}/{hap}", flush=True)
                continue
            arrays = junction_arrays(data) if spec.kind == "junctions" else {"values": np.asarray(data.values, dtype=np.float32)}
            resolution = None if spec.kind == "junctions" else int(getattr(data, "resolution", None) or spec.resolution or 1)
            if resolution is not None:
                arrays["resolution"] = np.asarray(resolution, dtype=np.int64)
            records = metadata_records(data)
            if not records:
                print(f"[WARN] {name} has no tracks for these tissues in {case.name}/{hap} (saved empty)", flush=True)
            npz_path = pred_dir / f"{spec.attr}.npz"
            tmp = pred_dir / f".{spec.attr}.tmp.npz"
            np.savez_compressed(tmp, **arrays)
            os.replace(tmp, npz_path)
            write_atomic_json(pred_dir / f"{spec.attr}_metadata.json", {"metadata": records, "output_type": name, "kind": spec.kind, "resolution": resolution})
            for record in records:
                curie = record.get("ontology_curie")
                if curie and curie not in self.details:
                    self.details[curie] = {k: record.get(k) for k in DETAIL_KEYS if k in record}
            (self.reference_outputs if hap == REFERENCE else self.done_outputs).setdefault(gene, set()).add(name)
            saved.append(spec.attr)
        return saved

    def update_metadata(self) -> None:
        if not self.done_outputs and not self.reference_outputs:
            return
        metadata = json.loads(self.meta_path.read_text(encoding="utf-8"))
        done = sorted({o for outputs in self.done_outputs.values() for o in outputs})
        if done:
            metadata["alphagenome_outputs"] = sorted(set(metadata.get("alphagenome_outputs") or []) | set(done))
        terms = list(self.ontology_terms or self.details.keys())
        metadata["ontologies"] = sorted(set(metadata.get("ontologies") or []) | set(terms))
        details = dict(metadata.get("ontology_details") or {})
        for curie, detail in self.details.items():
            details.setdefault(curie, detail)
        metadata["ontology_details"] = details
        catalog = metadata.get("window_catalog") or {}
        for gene in sorted(set(self.done_outputs) | set(self.reference_outputs)):
            window_path = self.dataset_dir / "references" / "windows" / gene / "window_metadata.json"
            entries = [catalog.setdefault(gene, {})]
            window_meta = None
            if window_path.exists():
                window_meta = json.loads(window_path.read_text(encoding="utf-8"))
                entries.append(window_meta)
            for entry in entries:
                if gene in self.done_outputs:
                    entry["outputs"] = sorted(set(entry.get("outputs") or []) | self.done_outputs[gene])
                if gene in self.reference_outputs:
                    entry["reference_outputs"] = sorted(set(entry.get("reference_outputs") or []) | self.reference_outputs[gene])
                entry["ontologies"] = sorted(set(entry.get("ontologies") or []) | set(terms))
            if window_meta is not None:
                write_atomic_json(window_path, window_meta)
        metadata["window_catalog"] = catalog
        metadata["last_updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        write_atomic_json(self.meta_path, metadata)

    def run(self) -> Dict[str, Any]:
        from alphagenome.models import dna_client

        from genomics.workflows.dataset_builders.non_longevous.build_window_and_predict import (
            _is_fatal_rpc_error,
            predict_sequence_resilient,
        )

        requested = {name: getattr(dna_client.OutputType, name) for name in self.outputs}
        todo = []
        skipped = 0
        if self.reference:
            for gene in self.genes:
                fasta, pending = self._pending(None, gene, REFERENCE)
                if fasta is None:
                    continue
                if pending:
                    todo.append(("reference", gene, REFERENCE, fasta, pending))
                else:
                    skipped += 1
        for sample in self.samples:
            for gene in self.genes:
                for hap in self.haplotypes:
                    fasta, pending = self._pending(sample, gene, hap)
                    if fasta is None:
                        continue
                    if pending:
                        todo.append((sample, gene, hap, fasta, pending))
                    else:
                        skipped += 1
        total = len(todo)
        print(f"[INFO] {total} haplotype windows to predict ({skipped} already have {', '.join(self.outputs)}); "
              f"ontologies: {', '.join(self.ontology_terms) if self.ontology_terms else 'ALL'}", flush=True)
        emit_progress(0.0, f"{total} haplotype windows to predict, {skipped} already done")
        if not total:
            print(f"@@result {json.dumps({'predicted': 0, 'skipped': skipped, 'failed': 0})}", flush=True)
            return {"predicted": 0, "skipped": skipped, "failed": 0}
        client_box = [self._client()]
        failures: List[str] = []
        consecutive = 0
        started = time.time()
        for index, (sample, gene, hap, fasta, pending) in enumerate(todo, start=1):
            sequence = read_fasta(fasta)
            if len(sequence) not in SUPPORTED_LENGTHS:
                failures.append(f"{sample}/{gene}/{hap}: length {len(sequence)} is not an AlphaGenome input length")
                print(f"[ERROR] {failures[-1]}", flush=True)
                continue
            try:
                prediction = predict_sequence_resilient(
                    client_box, self._client, sequence, [requested[o] for o in pending], self.ontology_terms,
                    timeout_s=self.timeout, max_attempts=self.max_attempts,
                )
                saved = self._save(fasta.parent, hap, prediction, pending, gene)
                consecutive = 0
                print(f"[INFO] {sample} {gene} {hap}: {', '.join(saved) or 'nothing saved'}", flush=True)
            except Exception as exc:
                if _is_fatal_rpc_error(exc):
                    self.update_metadata()
                    raise
                failures.append(f"{sample}/{gene}/{hap}: {type(exc).__name__}: {exc}")
                print(f"[ERROR] {failures[-1]}", flush=True)
                consecutive += 1
                if consecutive >= 5:
                    self.update_metadata()
                    raise RuntimeError(f"5 consecutive AlphaGenome failures; last: {exc}")
            elapsed = time.time() - started
            rate = elapsed / index
            eta = rate * (total - index)
            emit_progress(index / total, f"{sample} {gene} {hap} · {index}/{total} · {rate:.1f} s/window · ETA {eta / 60:.0f} min")
            if index % 25 == 0:
                self.update_metadata()
            if self.rate_limit_delay > 0:
                time.sleep(self.rate_limit_delay)
        self.update_metadata()
        result = {"predicted": total - len(failures), "skipped": skipped, "failed": len(failures)}
        print(f"@@result {json.dumps(result)}", flush=True)
        if failures and len(failures) == total:
            raise RuntimeError(f"Every prediction failed; first: {failures[0]}")
        return result


def _split(value: Optional[str]) -> List[str]:
    return [v.strip() for v in (value or "").split(",") if v.strip()]


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="genomics alphagenome predict-dataset", description="AlphaGenome predictions for every haplotype window of a canonical dataset")
    parser.add_argument("dataset_dir", type=Path)
    parser.add_argument("--outputs", default="RNA_SEQ", help=f"Comma-separated outputs among {', '.join(ALL_OUTPUTS)}, or 'all' (default: RNA_SEQ)")
    parser.add_argument("--ontology", default="", help="Comma-separated ontology CURIEs (tissues/cell types), e.g. CL:1000458,UBERON:0002107")
    parser.add_argument("--all-tissues", action="store_true", help="Allow an empty --ontology (every tissue: very large outputs)")
    parser.add_argument("--genes", default="", help="Comma-separated windows (default: all)")
    parser.add_argument("--samples", default="", help="Comma-separated sample ids (default: all)")
    parser.add_argument("--samples-file", type=Path, default=None, help="File with one sample id per line")
    parser.add_argument("--haplotypes", default="H1,H2", help="H1, H2 and/or ref (the reference window of each gene); 'ref' alone predicts only the reference")
    parser.add_argument("--overwrite", action="store_true", help="Re-predict windows that already have the outputs")
    parser.add_argument("--timeout", type=float, default=600.0, help="Per-call deadline in seconds (default: 600)")
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--rate-limit-delay", type=float, default=0.0, help="Seconds to wait after each call")
    args = parser.parse_args(argv)
    terms = _split(args.ontology)
    if not terms and not args.all_tissues:
        parser.error("give --ontology CURIEs (or --all-tissues to predict every tissue)")
    samples = _split(args.samples)
    if args.samples_file:
        samples += [line.strip() for line in args.samples_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    try:
        DatasetPredictor(
            args.dataset_dir,
            _split(args.outputs),
            terms,
            genes=_split(args.genes) or None,
            samples=samples or None,
            haplotypes=_split(args.haplotypes) or ["H1", "H2"],
            overwrite=args.overwrite,
            timeout=args.timeout,
            max_attempts=args.max_attempts,
            rate_limit_delay=args.rate_limit_delay,
        ).run()
    except (ValueError, FileNotFoundError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
