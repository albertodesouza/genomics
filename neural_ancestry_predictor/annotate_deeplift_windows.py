#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
O que o script faz:

1) Parseia seu TXT e extrai sequências FASTA (1000bp) por indivíduo/haplótipo/gene/center.
2) Baixa a sequência de referência hg38 (UCSC REST /getData/sequence).
   - A UCSC usa start 0-based (inclusive) e end 1-based (exclusive) para API/intervalos. :contentReference[oaicite:1]{index=1}
3) Chama variantes vs referência (SNVs por comparação; fallback com alinhamento simples se tamanhos divergirem).
4) Anota variantes com Ensembl VEP REST (POST /vep/homo_sapiens/region) em lotes (<=200). :contentReference[oaicite:2]{index=2}
5) Gera saídas agregadas:
   - por gene
   - por indivíduo (sample + hap)
   - por (indivíduo, gene)
6) Escreve um relatório Markdown com interpretação cautelosa para pigmentação.

Requisitos:
  pip install requests

Uso:
  python3 annotate_deeplift_windows.py <input.txt> --outdir out_variants

Observações:
- Este script depende de internet para UCSC e VEP REST.
- Para análises em larga escala/produção, considere rodar VEP local (script + cache), mas aqui mantemos REST.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable, Any

import requests


# ----------------------------
# Modelos
# ----------------------------

@dataclass
class RegionRecord:
    gene: str
    chrom: str
    center_1based: int
    sample: str
    hap: str
    header: str
    seq: str


@dataclass(frozen=True)
class VariantKey:
    chrom: str
    pos_1based: int
    ref: str
    alt: str


@dataclass
class VariantOcc:
    """Ocorrência de uma variante em um indivíduo/hap/gene/record."""
    key: VariantKey
    sample: str
    hap: str
    gene: str
    header: str
    center_1based: int  # Centro da janela (para calcular distância)
    distance_to_center: int  # Distância absoluta ao centro


@dataclass
class VariantAnnotation:
    key: VariantKey
    input_str: str
    most_severe_consequence: str
    worst_impact: str
    gene_symbols: str
    gene_ids: str
    transcript_ids: str
    protein_ids: str
    consequence_terms: str
    amino_acids: str
    codons: str
    rsids: str


# ----------------------------
# Regex de parsing (arquivo do usuário)
# ----------------------------

HEADER_RE = re.compile(
    r"^(?P<sample>[A-Za-z0-9]+)_(?P<hap>H[12])_(?P<gene>[A-Za-z0-9]+)_center_(?P<center>\d+)$"
)

ITEM_LINE_RE = re.compile(
    r"^\s*\d+\.\s*(?P<gene>[A-Za-z0-9]+)\s*:.*\b(?P<chrom>chr(?:[0-9]+|X|Y|M))\s*:\s*(?P<center>[\d,]+)\s*$"
)

DNA_BLOCK_RE = re.compile(
    r"^\s*DNA\s+(?P<hap>H[12])\s*\(1000bp\s+centradas\s+em\s+(?P<chrom>chr(?:[0-9]+|X|Y|M))\s*:\s*(?P<center>[\d,]+)\)\s*:\s*$",
    re.IGNORECASE
)


def parse_input_file(txt_path: Path) -> List[RegionRecord]:
    records: List[RegionRecord] = []
    current_gene: Optional[str] = None
    current_chrom: Optional[str] = None
    current_center: Optional[int] = None
    current_hap_from_block: Optional[str] = None

    lines = txt_path.read_text(encoding="utf-8", errors="replace").splitlines()

    i = 0
    while i < len(lines):
        ln = lines[i].strip()

        m_item = ITEM_LINE_RE.match(lines[i])
        if m_item:
            current_gene = m_item.group("gene")
            current_chrom = m_item.group("chrom")
            current_center = int(m_item.group("center").replace(",", ""))
            current_hap_from_block = None
            i += 1
            continue

        m_block = DNA_BLOCK_RE.match(lines[i])
        if m_block:
            current_hap_from_block = m_block.group("hap").upper()
            blk_chrom = m_block.group("chrom")
            blk_center = int(m_block.group("center").replace(",", ""))
            if current_chrom is None:
                current_chrom = blk_chrom
            if current_center is None:
                current_center = blk_center
            i += 1
            continue

        if ln.startswith(">"):
            header = ln[1:].strip()
            m_hdr = HEADER_RE.match(header)
            if not m_hdr:
                raise ValueError(f"Header FASTA inesperado: {header}")

            sample = m_hdr.group("sample")
            hap = m_hdr.group("hap")
            gene_from_hdr = m_hdr.group("gene")
            center_from_hdr = int(m_hdr.group("center"))

            gene = current_gene or gene_from_hdr
            chrom = current_chrom or "chr?"
            center = current_center or center_from_hdr

            if gene != gene_from_hdr:
                print(f"[WARN] Gene do item ({gene}) != gene do header ({gene_from_hdr}) em {header}", file=sys.stderr)
            if center != center_from_hdr:
                print(f"[WARN] Center do item ({center}) != center do header ({center_from_hdr}) em {header}", file=sys.stderr)
                center = center_from_hdr

            if current_hap_from_block and current_hap_from_block != hap:
                print(f"[WARN] Hap do bloco ({current_hap_from_block}) != hap do header ({hap}) em {header}", file=sys.stderr)

            seq_chunks: List[str] = []
            i += 1
            while i < len(lines):
                s = lines[i].strip()
                if not s:
                    break
                if s.startswith(">") or ITEM_LINE_RE.match(lines[i]) or DNA_BLOCK_RE.match(lines[i]):
                    break
                if re.fullmatch(r"[ACGTNacgtn]+", s):
                    seq_chunks.append(s.upper())
                i += 1

            seq = "".join(seq_chunks)
            if not seq:
                raise ValueError(f"Sequência vazia para header {header}")

            records.append(
                RegionRecord(
                    gene=gene,
                    chrom=chrom,
                    center_1based=center,
                    sample=sample,
                    hap=hap,
                    header=header,
                    seq=seq,
                )
            )
            continue

        i += 1

    if not records:
        raise RuntimeError("Nenhum record FASTA foi encontrado no arquivo de entrada.")
    return records


# ----------------------------
# UCSC referência hg38
# ----------------------------

def window_1000_centered(center_1based: int) -> Tuple[int, int]:
    """
    1000bp:
      center0 = center_1based - 1
      start0 = center0 - 500
      end1   = center0 + 500
    Comprimento = 1000 (end1 - start0)
    """
    center0 = center_1based - 1
    start0 = center0 - 500
    end1 = center0 + 500
    if start0 < 0:
        raise ValueError(f"Start < 0 para center={center_1based}")
    return start0, end1


def ucsc_get_sequence_hg38(chrom: str, start0: int, end1: int,
                          timeout_s: int = 30, retries: int = 5, backoff_s: float = 0.8) -> str:
    url = "https://api.genome.ucsc.edu/getData/sequence"
    params = {"genome": "hg38", "chrom": chrom, "start": str(start0), "end": str(end1)}
    headers = {"Accept": "application/json"}

    last_err: Optional[Exception] = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout_s)
            if r.status_code in (429, 503, 504):
                time.sleep(backoff_s * (2 ** attempt))
                continue
            r.raise_for_status()
            data = r.json()
            dna = data.get("dna")
            if not dna:
                raise RuntimeError(f"Resposta UCSC sem 'dna'. Keys: {list(data.keys())}")
            return str(dna).upper()
        except Exception as e:
            last_err = e
            time.sleep(backoff_s * (2 ** attempt))
    raise RuntimeError(f"Falha UCSC {chrom}:{start0}-{end1}. Erro: {last_err}")


def ref_cache_key(chrom: str, start0: int, end1: int) -> str:
    raw = f"{chrom}:{start0}-{end1}".encode("utf-8")
    return hashlib.sha1(raw).hexdigest()[:16]


# ----------------------------
# Chamada de variantes
# ----------------------------

def call_snvs_simple(rec: RegionRecord, ref_seq: str, start0: int) -> List[VariantKey]:
    if len(ref_seq) != len(rec.seq):
        raise ValueError("call_snvs_simple requer comprimentos iguais.")
    out: List[VariantKey] = []
    for i, (r, a) in enumerate(zip(ref_seq, rec.seq)):
        if r == a:
            continue
        pos1 = start0 + i + 1
        out.append(VariantKey(rec.chrom, pos1, r, a))
    return out


def needleman_wunsch(a: str, b: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> Tuple[str, str]:
    n, m = len(a), len(b)
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
    prev = [0] * (m + 1)
    cur = [0] * (m + 1)

    for j in range(1, m + 1):
        prev[j] = prev[j - 1] + gap
        ptr[0][j] = 2
    for i in range(1, n + 1):
        cur[0] = prev[0] + gap
        ptr[i][0] = 1
        for j in range(1, m + 1):
            s_diag = prev[j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            s_up = prev[j] + gap
            s_left = cur[j - 1] + gap
            best, best_ptr = s_diag, 0
            if s_up > best:
                best, best_ptr = s_up, 1
            if s_left > best:
                best, best_ptr = s_left, 2
            cur[j] = best
            ptr[i][j] = best_ptr
        prev, cur = cur, prev

    i, j = n, m
    out_a: List[str] = []
    out_b: List[str] = []
    while i > 0 or j > 0:
        p = ptr[i][j]
        if p == 0:
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif p == 1:
            out_a.append(a[i - 1]); out_b.append("-"); i -= 1
        else:
            out_a.append("-"); out_b.append(b[j - 1]); j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b))


def call_variants_with_alignment(rec: RegionRecord, ref_seq: str, start0: int) -> List[VariantKey]:
    """
    Fallback para indels; representação simples em key.ref/key.alt pode conter '-'.
    Se você quiser VCF estrito para indels, será necessário left-normalization com base âncora.
    """
    aln_ref, aln_alt = needleman_wunsch(ref_seq, rec.seq)
    out: List[VariantKey] = []

    ref_pos0 = start0
    i = 0
    while i < len(aln_ref):
        r = aln_ref[i]
        a = aln_alt[i]

        if r != "-" and a != "-":
            if r != a:
                out.append(VariantKey(rec.chrom, ref_pos0 + 1, r, a))
            ref_pos0 += 1
            i += 1
            continue

        if r != "-" and a == "-":  # deleção
            del_start0 = ref_pos0
            deleted = [r]
            ref_pos0 += 1
            i += 1
            while i < len(aln_ref) and aln_ref[i] != "-" and aln_alt[i] == "-":
                deleted.append(aln_ref[i])
                ref_pos0 += 1
                i += 1
            out.append(VariantKey(rec.chrom, del_start0 + 1, "".join(deleted), "-"))
            continue

        if r == "-" and a != "-":  # inserção
            ins_pos1 = ref_pos0  # aproximação
            inserted = [a]
            i += 1
            while i < len(aln_ref) and aln_ref[i] == "-" and aln_alt[i] != "-":
                inserted.append(aln_alt[i])
                i += 1
            out.append(VariantKey(rec.chrom, ins_pos1, "-", "".join(inserted)))
            continue

        i += 1

    return out


# ----------------------------
# Ensembl VEP REST (batching)
# ----------------------------

def chunked(lst: List[Any], n: int) -> Iterable[List[Any]]:
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def strip_chr(ch: str) -> str:
    return ch[3:] if ch.startswith("chr") else ch


def add_chr(ch: str) -> str:
    # Ensembl costuma devolver "11", "X", "Y", "MT"; padronizamos para UCSC "chr11", "chrX", ...
    if ch.startswith("chr"):
        return ch
    if ch in ("X", "Y", "M", "MT"):
        if ch == "MT":
            return "chrM"
        return f"chr{ch}"
    if ch.isdigit():
        return f"chr{ch}"
    return f"chr{ch}"


IMPACT_ORDER = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}


def worst_impact_from_transcripts(tcs: List[dict]) -> str:
    best = ""
    best_score = -1
    for tc in tcs or []:
        imp = str(tc.get("impact", "")).upper()
        score = IMPACT_ORDER.get(imp, 0)
        if score > best_score:
            best = imp
            best_score = score
    return best


def parse_vep_input_to_key(input_str: str) -> VariantKey:
    """
    Esperado (como enviamos): "11 61333312 . A G . . ."
    Pegamos tokens 0,1,3,4.
    """
    toks = input_str.strip().split()
    if len(toks) < 5:
        raise ValueError(f"Não consegui parsear VEP input: {input_str}")
    chrom = add_chr(toks[0])
    pos = int(toks[1])
    ref = toks[3].upper()
    alt = toks[4].upper()
    return VariantKey(chrom, pos, ref, alt)


def vep_annotate_variants_batched(
    unique_keys: List[VariantKey],
    outdir: Path,
    batch_size: int = 100,
    extra_params: Optional[Dict[str, str]] = None,
    timeout_s: int = 90,
    max_retries: int = 6,
    backoff_s: float = 1.0,
) -> List[dict]:
    """
    POST /vep/homo_sapiens/region
    Limite hard: 200 variantes por POST. :contentReference[oaicite:3]{index=3}
    """
    if batch_size > 200:
        raise ValueError("batch_size > 200 não permitido pelo endpoint VEP REST.")

    url = "https://rest.ensembl.org/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = extra_params or {}

    chunks_dir = outdir / "vep_chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)

    total = len(unique_keys)
    nbatches = ceil(total / batch_size)
    all_results: List[dict] = []

    for bi, batch in enumerate(chunked(unique_keys, batch_size), start=1):
        variant_lines = []
        for k in batch:
            # formato VCF-like por espaço: chrom pos . ref alt . . .
            variant_lines.append(f"{strip_chr(k.chrom)} {k.pos_1based} . {k.ref} {k.alt} . . .")
        payload = {"variants": variant_lines}

        last_err = None
        for attempt in range(max_retries):
            try:
                r = requests.post(url, headers=headers, params=params, data=json.dumps(payload), timeout=timeout_s)
                if r.status_code in (429, 503, 504):
                    time.sleep(backoff_s * (2 ** attempt))
                    continue
                if not r.ok:
                    raise RuntimeError(f"HTTP {r.status_code}: {r.text[:500]}")
                decoded = r.json()

                chunk_path = chunks_dir / f"vep_chunk_{bi:04d}_of_{nbatches:04d}.json"
                chunk_path.write_text(json.dumps(decoded, indent=2), encoding="utf-8")
                print(f"[SAVED] VEP chunk {bi}/{nbatches}: {chunk_path}")

                all_results.extend(decoded)
                break
            except Exception as e:
                last_err = e
                time.sleep(backoff_s * (2 ** attempt))
        else:
            raise RuntimeError(f"Falha definitiva no VEP (batch {bi}/{nbatches}). Último erro: {last_err}")

    aggregated_path = outdir / "vep_annotations.all.json"
    aggregated_path.write_text(json.dumps(all_results, indent=2), encoding="utf-8")
    print(f"[SAVED] VEP aggregated JSON: {aggregated_path}")

    return all_results


# ----------------------------
# TSV helpers
# ----------------------------

def write_tsv(path: Path, header: List[str], rows: List[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


# ----------------------------
# Agregação e relatório
# ----------------------------

PIGMENTATION_GENES = {
    "OCA2", "HERC2", "TYR", "SLC24A5", "SLC45A2", "MC1R", "MFSD12", "EDAR"
}

GENE_BACKGROUND = {
    "HERC2": "Lócus regulatório famoso que modula expressão de OCA2; variantes nessa região explicam grande parte da variação de cor dos olhos (e influenciam pigmentação). :contentReference[oaicite:4]{index=4}",
    "OCA2": "Gene fortemente associado a pigmentação; frequentemente co-analizado com HERC2 pela regulação enhancer-promoter (especialmente para cor dos olhos). :contentReference[oaicite:5]{index=5}",
    "SLC24A5": "Variante missense A111T (rs1426654) é um dos maiores determinantes de pele mais clara em populações europeias (grande efeito e alta frequência em europeus). :contentReference[oaicite:6]{index=6}",
    "SLC45A2": "Gene frequentemente associado à variação de cor de pele em humanos; aparece recorrentemente em modelos de pigmentação. :contentReference[oaicite:7]{index=7}",
    "TYR": "Enzima chave da melanogênese; variantes podem impactar produção de melanina (p.ex., quadros de hipopigmentação em alguns contextos).",
    "MC1R": "Receptor importante em via de melanogênese; variantes são clássicas em estudos de variação de pigmentação (p.ex., fenótipos de cabelo/pele em contextos populacionais).",
    "MFSD12": "Implicado em estudos genéticos de pigmentação; aparece como gene com associação a variações de cor de pele em análises populacionais.",
    "EDAR": "Lócus com variantes de grande efeito em traços ectodérmicos; frequentemente associado a assinaturas populacionais e, em alguns estudos, a traços relacionados."
}


def build_annotations(vep_json: List[dict]) -> Dict[VariantKey, VariantAnnotation]:
    out: Dict[VariantKey, VariantAnnotation] = {}

    for item in vep_json:
        inp = str(item.get("input", ""))
        if not inp:
            continue

        try:
            key = parse_vep_input_to_key(inp)
        except Exception:
            # fallback: tentar vcf_string se existir
            inp2 = str(item.get("vcf_string", "")) or inp
            key = parse_vep_input_to_key(inp2)

        msc = str(item.get("most_severe_consequence", ""))

        tcs = item.get("transcript_consequences") or []
        worst_imp = worst_impact_from_transcripts(tcs)

        gene_symbols = sorted({tc.get("gene_symbol","") for tc in tcs if tc.get("gene_symbol")})
        gene_ids = sorted({tc.get("gene_id","") for tc in tcs if tc.get("gene_id")})
        transcript_ids = sorted({tc.get("transcript_id","") for tc in tcs if tc.get("transcript_id")})
        protein_ids = sorted({tc.get("protein_id","") for tc in tcs if tc.get("protein_id")})

        # Consequence terms (agregando do primeiro TC, mas preservando um join útil)
        cons_terms = sorted({ct for tc in tcs for ct in (tc.get("consequence_terms") or [])})

        # Alguns campos convenientes do primeiro TC (se houver)
        amino_acids = ""
        codons = ""
        if tcs:
            amino_acids = str(tcs[0].get("amino_acids","") or "")
            codons = str(tcs[0].get("codons","") or "")

        # rsIDs (colocated_variants)
        rsids = sorted({cv.get("id","") for cv in (item.get("colocated_variants") or []) if cv.get("id")})

        out[key] = VariantAnnotation(
            key=key,
            input_str=inp,
            most_severe_consequence=msc,
            worst_impact=worst_imp,
            gene_symbols=",".join(gene_symbols),
            gene_ids=",".join(gene_ids),
            transcript_ids=",".join(transcript_ids),
            protein_ids=",".join(protein_ids),
            consequence_terms=",".join(cons_terms),
            amino_acids=amino_acids,
            codons=codons,
            rsids=",".join(rsids),
        )
    return out


def impact_bucket(imp: str) -> str:
    imp = (imp or "").upper()
    if imp in ("HIGH", "MODERATE", "LOW", "MODIFIER"):
        return imp
    return "UNKNOWN"


def summarize_gene_individual(
    occs: List[VariantOcc],
    ann: Dict[VariantKey, VariantAnnotation]
) -> Tuple[Dict[str, dict], Dict[str, dict], Dict[Tuple[str, str], dict]]:
    """
    Retorna:
      gene_summary[gene] = {...}
      indiv_summary[indiv] = {...} onde indiv = "SAMPLE_HAP"
      indiv_gene_summary[(indiv, gene)] = {...}
    """
    gene_sum: Dict[str, dict] = {}
    indiv_sum: Dict[str, dict] = {}
    indiv_gene_sum: Dict[Tuple[str, str], dict] = {}

    def bump(d: dict, k: str, v: int = 1) -> None:
        d[k] = int(d.get(k, 0)) + v

    def ensure(d: dict, k: str) -> dict:
        if k not in d:
            d[k] = {"total_occurrences": 0, "unique_variants": set(),
                    "impact_HIGH": 0, "impact_MODERATE": 0, "impact_LOW": 0, "impact_MODIFIER": 0, "impact_UNKNOWN": 0,
                    "top_msc": {}, "rsids": set()}
        return d[k]

    for o in occs:
        indiv = f"{o.sample}_{o.hap}"
        a = ann.get(o.key)

        gene_entry = ensure(gene_sum, o.gene)
        indiv_entry = ensure(indiv_sum, indiv)
        ig_entry = ensure(indiv_gene_sum, (indiv, o.gene))

        for entry in (gene_entry, indiv_entry, ig_entry):
            entry["total_occurrences"] += 1
            entry["unique_variants"].add(o.key)

        if a:
            b = impact_bucket(a.worst_impact)
            for entry in (gene_entry, indiv_entry, ig_entry):
                bump(entry, f"impact_{b}", 1)
                if a.most_severe_consequence:
                    entry["top_msc"][a.most_severe_consequence] = entry["top_msc"].get(a.most_severe_consequence, 0) + 1
                if a.rsids:
                    for rid in a.rsids.split(","):
                        if rid:
                            entry["rsids"].add(rid)
        else:
            for entry in (gene_entry, indiv_entry, ig_entry):
                bump(entry, "impact_UNKNOWN", 1)

    # converter sets -> contagens
    for d in list(gene_sum.values()) + list(indiv_sum.values()) + list(indiv_gene_sum.values()):
        d["unique_variants_count"] = len(d["unique_variants"])
        del d["unique_variants"]

        # top msc
        top_msc = sorted(d["top_msc"].items(), key=lambda x: (-x[1], x[0]))[:5]
        d["top_msc_str"] = ",".join([f"{k}({v})" for k, v in top_msc])
        del d["top_msc"]

        # rsids
        d["rsids_str"] = ",".join(sorted(d["rsids"])) if d["rsids"] else ""
        del d["rsids"]

    return gene_sum, indiv_sum, indiv_gene_sum


def write_report_markdown(
    report_path: Path,
    gene_sum: Dict[str, dict],
    indiv_sum: Dict[str, dict],
    indiv_gene_sum: Dict[Tuple[str, str], dict],
    ann: Dict[VariantKey, VariantAnnotation],
    occs: List[VariantOcc]
) -> None:
    """
    Gera um relatório textual com:
      - visão geral
      - destaques por gene
      - destaques por indivíduo
      - interpretação cautelosa para pigmentação (com base em genes clássicos)
    """
    n_occ = len(occs)
    uniq_vars = len({o.key for o in occs})

    genes = sorted(gene_sum.keys())
    indivs = sorted(indiv_sum.keys())

    # top genes por HIGH+MODERATE (proxy de efeitos em proteína/splicing/regulação anotada)
    def score_gene(g: str) -> Tuple[int, int, int]:
        d = gene_sum[g]
        return (d.get("impact_HIGH", 0) + d.get("impact_MODERATE", 0),
                d.get("impact_HIGH", 0),
                d.get("unique_variants_count", 0))

    top_genes = sorted(genes, key=score_gene, reverse=True)[:15]

    # top indiv por HIGH+MODERATE
    def score_ind(i: str) -> Tuple[int, int, int]:
        d = indiv_sum[i]
        return (d.get("impact_HIGH", 0) + d.get("impact_MODERATE", 0),
                d.get("impact_HIGH", 0),
                d.get("unique_variants_count", 0))

    top_indivs = sorted(indivs, key=score_ind, reverse=True)[:15]

    # helper: gene narrative
    def gene_context(g: str) -> str:
        if g in GENE_BACKGROUND:
            return GENE_BACKGROUND[g]
        return ""

    # pigmentação: explicação geral + cautelas
    pigment_intro = (
        "### Interpretação cautelosa: o que dá (e o que não dá) para inferir\n"
        "- Pigmentação (pele/olhos/cabelo) é **poligênica** e altamente dependente do contexto populacional; variantes em poucos lócus podem ter grande efeito "
        "em certas populações (ex.: HERC2/OCA2 para cor dos olhos; SLC24A5 em variação de pele mais clara em europeus). :contentReference[oaicite:8]{index=8}\n"
        "- O VEP indica **consequências funcionais** (p.ex., missense, splice, intronic) e frequentemente retorna **rsIDs** (quando a variante coincide com uma já catalogada). "
        "Isso ajuda a priorizar candidatos, mas **não prova** aumento/diminuição de expressão ou “produção maior de proteína” sem dados adicionais (eQTL/expressão/proteômica).\n"
        "- O que este relatório faz é: **destacar** genes e indivíduos com (i) maior carga de variantes, (ii) consequências mais severas e (iii) presença de rsIDs conhecidos, "
        "especialmente em genes clássicos de pigmentação.\n"
    )

    with report_path.open("w", encoding="utf-8") as f:
        f.write("# Variant Annotation Report (por gene / por indivíduo)\n\n")
        f.write("## Visão geral\n")
        f.write(f"- Ocorrências (indivíduo×janela): **{n_occ}**\n")
        f.write(f"- Variantes únicas (chrom,pos,ref,alt): **{uniq_vars}**\n")
        f.write(f"- Genes (do seu arquivo): **{len(genes)}**\n")
        f.write(f"- Indivíduos/hap (SAMPLE_HAP): **{len(indivs)}**\n\n")

        f.write("## Destaques por gene (ordenado por HIGH+MODERATE)\n\n")
        f.write("| Gene | Unique variants | Occurrences | HIGH | MODERATE | LOW | MODIFIER | Top consequences | rsIDs |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---|---|\n")
        for g in top_genes:
            d = gene_sum[g]
            f.write(
                f"| {g} | {d['unique_variants_count']} | {d['total_occurrences']} | "
                f"{d.get('impact_HIGH',0)} | {d.get('impact_MODERATE',0)} | {d.get('impact_LOW',0)} | {d.get('impact_MODIFIER',0)} | "
                f"{d.get('top_msc_str','')} | {d.get('rsids_str','')} |\n"
            )
        f.write("\n")

        # contexto para genes de pigmentação, se presentes
        present_pig = [g for g in genes if g in PIGMENTATION_GENES]
        if present_pig:
            f.write("### Genes clássicos de pigmentação presentes no seu conjunto\n")
            for g in present_pig:
                d = gene_sum[g]
                f.write(f"- **{g}**: unique={d['unique_variants_count']}, HIGH={d.get('impact_HIGH',0)}, MODERATE={d.get('impact_MODERATE',0)}. ")
                ctx = gene_context(g)
                if ctx:
                    f.write(ctx)
                f.write("\n")
            f.write("\n")
        else:
            f.write("### Genes clássicos de pigmentação\n")
            f.write("Nenhum dos genes clássicos (OCA2/HERC2/SLC24A5/SLC45A2/TYR/MC1R/MFSD12/EDAR) apareceu explicitamente no seu arquivo. "
                    "Se você esperava esses genes, verifique a lista de genes/janelas de entrada.\n\n")

        f.write("## Destaques por indivíduo (SAMPLE_HAP) (ordenado por HIGH+MODERATE)\n\n")
        f.write("| Indivíduo | Unique variants | Occurrences | HIGH | MODERATE | LOW | MODIFIER | Top consequences | rsIDs |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---|---|\n")
        for i_id in top_indivs:
            d = indiv_sum[i_id]
            f.write(
                f"| {i_id} | {d['unique_variants_count']} | {d['total_occurrences']} | "
                f"{d.get('impact_HIGH',0)} | {d.get('impact_MODERATE',0)} | {d.get('impact_LOW',0)} | {d.get('impact_MODIFIER',0)} | "
                f"{d.get('top_msc_str','')} | {d.get('rsids_str','')} |\n"
            )
        f.write("\n")

        f.write("## Matriz indivíduo × gene (resumo)\n\n")
        f.write("Abaixo estão os pares (indivíduo, gene) com maior carga de HIGH+MODERATE.\n\n")
        pairs = list(indiv_gene_sum.keys())

        def score_pair(p: Tuple[str, str]) -> Tuple[int, int, int]:
            d = indiv_gene_sum[p]
            return (d.get("impact_HIGH", 0) + d.get("impact_MODERATE", 0),
                    d.get("impact_HIGH", 0),
                    d.get("unique_variants_count", 0))

        top_pairs = sorted(pairs, key=score_pair, reverse=True)[:30]
        f.write("| Indivíduo | Gene | Unique variants | Occurrences | HIGH | MODERATE | Top consequences | rsIDs |\n")
        f.write("|---|---|---:|---:|---:|---:|---|---|\n")
        for (i_id, g) in top_pairs:
            d = indiv_gene_sum[(i_id, g)]
            f.write(
                f"| {i_id} | {g} | {d['unique_variants_count']} | {d['total_occurrences']} | "
                f"{d.get('impact_HIGH',0)} | {d.get('impact_MODERATE',0)} | {d.get('top_msc_str','')} | {d.get('rsids_str','')} |\n"
            )
        f.write("\n")

        f.write(pigment_intro)

        # Inferência por indivíduo focada em genes clássicos (se existirem)
        if present_pig:
            f.write("### Interpretação por indivíduo (focada em genes de pigmentação)\n")
            f.write(
                "Para cada indivíduo, observe especialmente:\n"
                "- presença de **rsIDs** bem estudados (quando aparecem);\n"
                "- consequências **MODERATE/HIGH** em genes de pigmentação;\n"
                "- padrões consistentes entre H1 e H2.\n\n"
            )
            for i_id in sorted(indivs)[:50]:
                # resumir somente se indivíduo tocar genes de pigmentação
                touched = [(g, indiv_gene_sum[(i_id, g)]) for g in present_pig if (i_id, g) in indiv_gene_sum]
                if not touched:
                    continue
                f.write(f"#### {i_id}\n")
                for g, d in sorted(touched, key=lambda x: (-(x[1].get("impact_HIGH",0)+x[1].get("impact_MODERATE",0)), x[0])):
                    f.write(
                        f"- {g}: unique={d['unique_variants_count']}, HIGH={d.get('impact_HIGH',0)}, MODERATE={d.get('impact_MODERATE',0)}, "
                        f"top={d.get('top_msc_str','')}, rsIDs={d.get('rsids_str','')}\n"
                    )
                f.write(
                    "Interpretação: se houver rsIDs clássicos e/ou consequências MODERATE/HIGH em genes como HERC2/OCA2 (cor dos olhos) ou SLC24A5/SLC45A2 (pele), "
                    "isso pode indicar um perfil genético com maior probabilidade de certas tendências de pigmentação, mas a conclusão depende de quais alelos específicos "
                    "foram observados e do contexto populacional. :contentReference[oaicite:9]{index=9}\n\n"
                )

        # Nota final com pointers
        f.write("## Notas práticas para aprofundar\n")
        f.write(
            "1) Se você quer inferir “produção maior de proteína”, priorize variantes em regiões **codificantes** (missense/nonsense/frameshift) e em **splicing**; "
            "para efeitos de expressão, cruze com eQTL/regulatory tracks.\n"
            "2) Se o VEP retornar rsIDs (colocated_variants), você pode cruzar com literatura/GTEx/ClinVar/Ensembl Variation.\n"
            "3) Para lotes muito grandes, VEP local (script + cache) tende a ser mais estável do que o REST. :contentReference[oaicite:10]{index=10}\n"
        )


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input_txt", type=str, help="Arquivo TXT do seu pipeline (com FASTAs).")
    ap.add_argument("--outdir", type=str, default="variant_annotation_out", help="Diretório de saída.")
    ap.add_argument("--vep-batch-size", type=int, default=100, help="Batch VEP (<=200).")
    ap.add_argument("--ucsc-timeout", type=int, default=30)
    ap.add_argument("--vep-timeout", type=int, default=90)
    ap.add_argument("--no-vep", action="store_true", help="Pula VEP (gera só variantes vs referência).")
    ap.add_argument(
        "--central-window", 
        type=int, 
        default=None, 
        metavar="BP",
        help="Filtrar variantes para apenas aquelas dentro de ±(BP/2) do centro. "
             "Ex: --central-window 50 mantém apenas variantes dentro de ±25bp do centro."
    )
    args = ap.parse_args()

    input_path = Path(args.input_txt).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    refs_dir = outdir / "reference_fasta"
    refs_dir.mkdir(exist_ok=True)

    # Outputs principais
    parsed_records_path = outdir / "parsed_records.tsv"
    occs_path = outdir / "variant_occurrences.tsv"
    uniq_vars_path = outdir / "unique_variants.tsv"
    vep_ann_path = outdir / "vep_annotated_variants.tsv"
    gene_sum_path = outdir / "summary_by_gene.tsv"
    indiv_sum_path = outdir / "summary_by_individual.tsv"
    indiv_gene_sum_path = outdir / "summary_by_individual_gene.tsv"
    report_path = outdir / "report_by_gene_by_individual.md"

    records = parse_input_file(input_path)

    # Salva parsed records
    pr_header = ["gene", "chrom", "center_1based", "sample", "hap", "header", "seq_len"]
    pr_rows = [[r.gene, r.chrom, str(r.center_1based), r.sample, r.hap, r.header, str(len(r.seq))] for r in records]
    write_tsv(parsed_records_path, pr_header, pr_rows)
    print(f"[SAVED] parsed records TSV: {parsed_records_path}")

    # Cache de referência por região
    ref_cache: Dict[str, str] = {}

    # Coletar ocorrências e variantes únicas
    occs: List[VariantOcc] = []
    unique_keys_set: Dict[VariantKey, None] = {}

    for rec in records:
        start0, end1 = window_1000_centered(rec.center_1based)
        cache_k = ref_cache_key(rec.chrom, start0, end1)

        if cache_k in ref_cache:
            ref_seq = ref_cache[cache_k]
        else:
            ref_seq = ucsc_get_sequence_hg38(
                rec.chrom, start0, end1,
                timeout_s=args.ucsc_timeout,
                retries=5, backoff_s=0.8
            )
            ref_cache[cache_k] = ref_seq

            ref_fa = refs_dir / f"{rec.chrom}_{start0}_{end1}.ref.fa"
            if not ref_fa.exists():
                with ref_fa.open("w", encoding="utf-8") as f:
                    f.write(f">{rec.chrom}:{start0}-{end1}(hg38)\n")
                    for k in range(0, len(ref_seq), 60):
                        f.write(ref_seq[k:k+60] + "\n")

        # Chamada de variantes
        if len(ref_seq) == len(rec.seq):
            keys = call_snvs_simple(rec, ref_seq, start0)
        else:
            keys = call_variants_with_alignment(rec, ref_seq, start0)

        for k in keys:
            occs.append(VariantOcc(key=k, sample=rec.sample, hap=rec.hap, gene=rec.gene, header=rec.header))
            unique_keys_set[k] = None

    if not occs:
        print("[INFO] Nenhuma variante detectada. Encerrando.")
        return 0

    unique_keys = list(unique_keys_set.keys())
    print(f"[INFO] Variant occurrences: {len(occs)} | unique variants: {len(unique_keys)}")

    # Salvar ocorrências
    occ_header = ["chrom", "pos_1based", "ref", "alt", "sample", "hap", "gene", "header"]
    occ_rows = [[o.key.chrom, str(o.key.pos_1based), o.key.ref, o.key.alt, o.sample, o.hap, o.gene, o.header] for o in occs]
    write_tsv(occs_path, occ_header, occ_rows)
    print(f"[SAVED] variant occurrences TSV: {occs_path}")

    # Salvar variantes únicas
    uv_header = ["chrom", "pos_1based", "ref", "alt"]
    uv_rows = [[k.chrom, str(k.pos_1based), k.ref, k.alt] for k in sorted(unique_keys, key=lambda x: (x.chrom, x.pos_1based, x.ref, x.alt))]
    write_tsv(uniq_vars_path, uv_header, uv_rows)
    print(f"[SAVED] unique variants TSV: {uniq_vars_path}")

    # Rodar VEP (se habilitado)
    ann_map: Dict[VariantKey, VariantAnnotation] = {}
    if not args.no_vep:
        vep_params = {
            "hgvs": "1",
            "protein": "1",
            "uniprot": "1",
            "canonical": "1",
            "mane": "1",
            "vcf_string": "1",
        }

        vep_json = vep_annotate_variants_batched(
            unique_keys,
            outdir=outdir,
            batch_size=args.vep_batch_size,
            extra_params=vep_params,
            timeout_s=args.vep_timeout,
            max_retries=6,
            backoff_s=1.0
        )

        ann_map = build_annotations(vep_json)

        # Salvar anotações (TSV)
        va_header = [
            "chrom", "pos_1based", "ref", "alt",
            "most_severe_consequence", "worst_impact",
            "gene_symbols", "gene_ids", "transcript_ids", "protein_ids",
            "consequence_terms", "amino_acids", "codons", "rsids"
        ]
        va_rows: List[List[str]] = []
        for k, a in sorted(ann_map.items(), key=lambda x: (x[0].chrom, x[0].pos_1based, x[0].ref, x[0].alt)):
            va_rows.append([
                k.chrom, str(k.pos_1based), k.ref, k.alt,
                a.most_severe_consequence, a.worst_impact,
                a.gene_symbols, a.gene_ids, a.transcript_ids, a.protein_ids,
                a.consequence_terms, a.amino_acids, a.codons, a.rsids
            ])
        write_tsv(vep_ann_path, va_header, va_rows)
        print(f"[SAVED] VEP annotated variants TSV: {vep_ann_path}")
    else:
        print("[INFO] --no-vep habilitado: pulando anotação VEP.")

    # Agregações (usam ann_map quando existe; caso não, buckets UNKNOWN)
    gene_sum, indiv_sum, indiv_gene_sum = summarize_gene_individual(occs, ann_map)

    # write summary_by_gene
    g_header = ["gene", "unique_variants", "occurrences", "HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN", "top_consequences", "rsids"]
    g_rows: List[List[str]] = []
    for g, d in sorted(gene_sum.items(), key=lambda x: (-(x[1].get("impact_HIGH",0)+x[1].get("impact_MODERATE",0)), x[0])):
        g_rows.append([
            g,
            str(d["unique_variants_count"]),
            str(d["total_occurrences"]),
            str(d.get("impact_HIGH",0)),
            str(d.get("impact_MODERATE",0)),
            str(d.get("impact_LOW",0)),
            str(d.get("impact_MODIFIER",0)),
            str(d.get("impact_UNKNOWN",0)),
            d.get("top_msc_str",""),
            d.get("rsids_str","")
        ])
    write_tsv(gene_sum_path, g_header, g_rows)
    print(f"[SAVED] summary by gene TSV: {gene_sum_path}")

    # write summary_by_individual
    i_header = ["individual", "unique_variants", "occurrences", "HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN", "top_consequences", "rsids"]
    i_rows: List[List[str]] = []
    for ind, d in sorted(indiv_sum.items(), key=lambda x: (-(x[1].get("impact_HIGH",0)+x[1].get("impact_MODERATE",0)), x[0])):
        i_rows.append([
            ind,
            str(d["unique_variants_count"]),
            str(d["total_occurrences"]),
            str(d.get("impact_HIGH",0)),
            str(d.get("impact_MODERATE",0)),
            str(d.get("impact_LOW",0)),
            str(d.get("impact_MODIFIER",0)),
            str(d.get("impact_UNKNOWN",0)),
            d.get("top_msc_str",""),
            d.get("rsids_str","")
        ])
    write_tsv(indiv_sum_path, i_header, i_rows)
    print(f"[SAVED] summary by individual TSV: {indiv_sum_path}")

    # write summary_by_individual_gene
    ig_header = ["individual", "gene", "unique_variants", "occurrences", "HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN", "top_consequences", "rsids"]
    ig_rows: List[List[str]] = []
    for (ind, g), d in sorted(indiv_gene_sum.items(), key=lambda x: (-(x[1].get("impact_HIGH",0)+x[1].get("impact_MODERATE",0)), x[0][0], x[0][1])):
        ig_rows.append([
            ind, g,
            str(d["unique_variants_count"]),
            str(d["total_occurrences"]),
            str(d.get("impact_HIGH",0)),
            str(d.get("impact_MODERATE",0)),
            str(d.get("impact_LOW",0)),
            str(d.get("impact_MODIFIER",0)),
            str(d.get("impact_UNKNOWN",0)),
            d.get("top_msc_str",""),
            d.get("rsids_str","")
        ])
    write_tsv(indiv_gene_sum_path, ig_header, ig_rows)
    print(f"[SAVED] summary by individual+gene TSV: {indiv_gene_sum_path}")

    # relatório textual
    write_report_markdown(report_path, gene_sum, indiv_sum, indiv_gene_sum, ann_map, occs)
    print(f"[SAVED] report (Markdown): {report_path}")

    print("\n=== OUTPUTS PRINCIPAIS ===")
    print(f"[OUT] {parsed_records_path}")
    print(f"[OUT] {refs_dir} (FASTA de referência por região)")
    print(f"[OUT] {occs_path}")
    print(f"[OUT] {uniq_vars_path}")
    if not args.no_vep:
        print(f"[OUT] {outdir / 'vep_chunks'} (JSON por lote)")
        print(f"[OUT] {outdir / 'vep_annotations.all.json'}")
        print(f"[OUT] {vep_ann_path}")
    print(f"[OUT] {gene_sum_path}")
    print(f"[OUT] {indiv_sum_path}")
    print(f"[OUT] {indiv_gene_sum_path}")
    print(f"[OUT] {report_path}")
    print("==========================\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

