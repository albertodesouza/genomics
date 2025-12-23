#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
O que o script faz:

1) Parseia seu TXT e extrai sequências FASTA por indivíduo/haplótipo/gene/center.
   - O tamanho da janela é detectado automaticamente do comprimento das sequências.
2) Baixa a sequência de referência hg38 (UCSC REST /getData/sequence).
   - A UCSC usa start 0-based (inclusive) e end 1-based (exclusive) para API/intervalos.
3) Chama variantes vs referência (SNVs por comparação; fallback com alinhamento simples se tamanhos divergirem).
4) Anota variantes com Ensembl VEP REST (POST /vep/homo_sapiens/region) em lotes (<=200).
5) Gera saídas agregadas:
   - por gene
   - por indivíduo (sample + hap)
   - por (indivíduo, gene)
6) Escreve um relatório Markdown genérico para validação de fenótipo.

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
    r"^\s*DNA\s+(?P<hap>H[12])\s*\((?P<bp>\d+)bp\s+centradas\s+em\s+(?P<chrom>chr(?:[0-9]+|X|Y|M))\s*:\s*(?P<center>[\d,]+)\)\s*:\s*$",
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

def window_centered(center_1based: int, window_size: int) -> Tuple[int, int]:
    """
    Retorna (start0, end1) para janela de window_size bp centrada em center_1based.
    
    Args:
        center_1based: Posição central (1-based)
        window_size: Tamanho da janela em bp
    
    Returns:
        Tupla (start0, end1) onde start0 é 0-based e end1 é exclusivo
    """
    center0 = center_1based - 1
    half = window_size // 2
    start0 = center0 - half
    end1 = center0 + (window_size - half)  # Garante tamanho exato para janelas ímpares
    if start0 < 0:
        raise ValueError(f"Start < 0 para center={center_1based}, window_size={window_size}")
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
# Gene Information APIs
# ----------------------------

@dataclass
class GeneInfo:
    """Information about a gene fetched from external APIs."""
    symbol: str
    description: str
    function: str
    biotype: str
    chromosome: str
    strand: str  # "+" or "-" for protein expression strand
    source: str


def _fetch_ensembl_gene(gene_symbol: str, timeout_s: int = 15) -> Optional[dict]:
    """Fetch gene info from Ensembl REST API."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    try:
        r = requests.get(url, headers=headers, timeout=timeout_s)
        if r.status_code == 200:
            return r.json()
    except Exception:
        pass
    return None


def _fetch_ncbi_gene(gene_symbol: str, timeout_s: int = 15) -> Optional[dict]:
    """Fetch gene info from NCBI Gene API."""
    # First, search for the gene
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "gene",
        "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
        "retmode": "json",
        "retmax": "1"
    }
    
    try:
        r = requests.get(search_url, params=search_params, timeout=timeout_s)
        if r.status_code != 200:
            return None
        
        data = r.json()
        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None
        
        gene_id = id_list[0]
        
        # Now fetch summary
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "json"
        }
        
        r = requests.get(summary_url, params=summary_params, timeout=timeout_s)
        if r.status_code == 200:
            data = r.json()
            result = data.get("result", {})
            if gene_id in result:
                return result[gene_id]
    except Exception:
        pass
    return None


def _fetch_uniprot_gene(gene_symbol: str, timeout_s: int = 15) -> Optional[dict]:
    """Fetch protein function info from UniProt API."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": "1",
        "fields": "gene_names,protein_name,cc_function"
    }
    
    try:
        r = requests.get(url, params=params, timeout=timeout_s)
        if r.status_code == 200:
            data = r.json()
            results = data.get("results", [])
            if results:
                return results[0]
    except Exception:
        pass
    return None


def fetch_gene_info(
    gene_symbols: List[str], 
    rate_limit_s: float = 0.34,  # ~3 req/s to stay under limits
    timeout_s: int = 15
) -> Dict[str, GeneInfo]:
    """
    Fetch information about genes from multiple APIs.
    
    Combines data from:
    - Ensembl: description, biotype, chromosome
    - NCBI Gene: functional summary
    - UniProt: protein function
    
    Returns a dictionary mapping gene symbols to GeneInfo objects.
    """
    gene_info_map: Dict[str, GeneInfo] = {}
    
    for gene in gene_symbols:
        description = ""
        function = ""
        biotype = ""
        chromosome = ""
        strand = ""
        sources = []
        
        # Ensembl
        ensembl_data = _fetch_ensembl_gene(gene, timeout_s)
        if ensembl_data:
            description = ensembl_data.get("description", "") or ""
            biotype = ensembl_data.get("biotype", "") or ""
            chromosome = ensembl_data.get("seq_region_name", "") or ""
            # Strand: 1 = forward (+), -1 = reverse (-)
            strand_int = ensembl_data.get("strand")
            if strand_int == 1:
                strand = "+"
            elif strand_int == -1:
                strand = "-"
            sources.append("Ensembl")
        
        time.sleep(rate_limit_s)
        
        # NCBI
        ncbi_data = _fetch_ncbi_gene(gene, timeout_s)
        if ncbi_data:
            ncbi_desc = ncbi_data.get("description", "") or ""
            ncbi_summary = ncbi_data.get("summary", "") or ""
            if ncbi_summary:
                # Truncate long summaries
                if len(ncbi_summary) > 300:
                    ncbi_summary = ncbi_summary[:297] + "..."
                function = ncbi_summary
            if not description and ncbi_desc:
                description = ncbi_desc
            sources.append("NCBI")
        
        time.sleep(rate_limit_s)
        
        # UniProt
        uniprot_data = _fetch_uniprot_gene(gene, timeout_s)
        if uniprot_data:
            # Extract function from comments
            comments = uniprot_data.get("comments", [])
            for comment in comments:
                if comment.get("commentType") == "FUNCTION":
                    texts = comment.get("texts", [])
                    if texts:
                        uniprot_func = texts[0].get("value", "")
                        if uniprot_func and (not function or len(uniprot_func) > len(function)):
                            if len(uniprot_func) > 300:
                                uniprot_func = uniprot_func[:297] + "..."
                            function = uniprot_func
            sources.append("UniProt")
        
        time.sleep(rate_limit_s)
        
        # Create GeneInfo
        gene_info_map[gene] = GeneInfo(
            symbol=gene,
            description=description or f"Gene {gene}",
            function=function or "Function not available from APIs.",
            biotype=biotype or "unknown",
            chromosome=chromosome or "unknown",
            strand=strand or "unknown",
            source=", ".join(sources) if sources else "No data found"
        )
        
        strand_str = f", strand={strand}" if strand else ""
        print(f"[API] Gene {gene}: fetched from {', '.join(sources) if sources else 'no sources'}{strand_str}")
    
    return gene_info_map


# ----------------------------
# rsID Information API
# ----------------------------

@dataclass
class RsidInfo:
    """Information about an rsID fetched from Ensembl Variation."""
    rsid: str
    description: str
    clinical_significance: str
    phenotypes: str
    minor_allele: str
    maf: str
    source: str


# ----------------------------
# HIGH Impact Variant Details
# ----------------------------

@dataclass
class HighImpactVariantInfo:
    """Detailed information about a HIGH impact variant."""
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    consequence: str
    amino_acids: str  # e.g., "Y/*" for stop_gained
    codons: str  # e.g., "taC/taG"
    protein_position: str
    transcript_id: str
    rsid: str
    clinical_significance: str
    phenotypes: str
    cosmic_id: str
    maf: str
    distance_to_center: int
    ensembl_link: str
    source: str


def fetch_high_impact_details(
    high_variants: List[Tuple['VariantKey', 'VariantAnnotation', int]],
    rate_limit_s: float = 0.34,
    timeout_s: int = 15
) -> List[HighImpactVariantInfo]:
    """
    Fetch detailed information about HIGH impact variants from Ensembl APIs.
    
    Args:
        high_variants: List of (VariantKey, VariantAnnotation, distance_to_center) tuples
        rate_limit_s: Rate limit between API calls
        timeout_s: Timeout for API calls
    
    Returns:
        List of HighImpactVariantInfo objects with detailed information
    """
    results: List[HighImpactVariantInfo] = []
    
    for key, ann, distance in high_variants:
        # Extract basic info from annotation
        gene = ann.gene_symbols.split(",")[0] if ann.gene_symbols else "Unknown"
        consequence = ann.most_severe_consequence
        amino_acids = ann.amino_acids or ""
        codons = ann.codons or ""
        transcript_id = ann.transcript_ids.split(",")[0] if ann.transcript_ids else ""
        
        # Get rsID if available
        rsid = ""
        if ann.rsids:
            rsids = [r for r in ann.rsids.split(",") if r.startswith("rs")]
            if rsids:
                rsid = rsids[0]
        
        # Get COSMIC ID if available
        cosmic_id = ""
        if ann.rsids:
            cosmics = [r for r in ann.rsids.split(",") if r.startswith("COSV")]
            if cosmics:
                cosmic_id = cosmics[0]
        
        clinical_sig = ""
        phenotypes_str = ""
        maf = ""
        protein_position = ""
        
        # Try to get more info from Ensembl Variation if we have an rsID
        if rsid:
            url = f"https://rest.ensembl.org/variation/human/{rsid}"
            headers = {"Content-Type": "application/json", "Accept": "application/json"}
            
            try:
                r = requests.get(url, headers=headers, timeout=timeout_s)
                if r.status_code == 200:
                    data = r.json()
                    
                    # Clinical significance
                    clinical_sigs = data.get("clinical_significance", [])
                    if clinical_sigs:
                        clinical_sig = ", ".join(clinical_sigs)
                    
                    # Phenotypes
                    phenotypes = data.get("phenotypes", [])
                    if phenotypes:
                        pheno_names = [p.get("trait", "") for p in phenotypes if p.get("trait")]
                        if pheno_names:
                            phenotypes_str = "; ".join(pheno_names[:3])  # Top 3
                            if len(pheno_names) > 3:
                                phenotypes_str += f" (+{len(pheno_names)-3} more)"
                    
                    # MAF
                    maf_val = data.get("MAF")
                    if maf_val is not None:
                        maf = f"{maf_val:.4f}"
                    
                    print(f"[API] HIGH variant {rsid}: clinical={clinical_sig or 'N/A'}")
                    
            except Exception as e:
                print(f"[WARN] Failed to fetch details for {rsid}: {e}", file=sys.stderr)
            
            time.sleep(rate_limit_s)
        
        # Try to get protein position from VEP consequence terms
        # Parse from annotation if available (format: "position/total")
        if ann.consequence_terms:
            # VEP sometimes includes position info
            pass
        
        # Build Ensembl link
        ensembl_link = f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={rsid}" if rsid else ""
        if not ensembl_link and cosmic_id:
            ensembl_link = f"https://cancer.sanger.ac.uk/cosmic/search?q={cosmic_id}"
        if not ensembl_link:
            # Link to region
            ensembl_link = f"https://www.ensembl.org/Homo_sapiens/Location/View?r={key.chrom.replace('chr', '')}:{key.pos_1based-50}-{key.pos_1based+50}"
        
        results.append(HighImpactVariantInfo(
            chrom=key.chrom,
            pos=key.pos_1based,
            ref=key.ref,
            alt=key.alt,
            gene=gene,
            consequence=consequence,
            amino_acids=amino_acids,
            codons=codons,
            protein_position=protein_position,
            transcript_id=transcript_id,
            rsid=rsid,
            clinical_significance=clinical_sig or "Not reported",
            phenotypes=phenotypes_str or "No phenotypes reported",
            cosmic_id=cosmic_id,
            maf=maf or "N/A",
            distance_to_center=distance,
            ensembl_link=ensembl_link,
            source="Ensembl Variation" if rsid else "VEP annotation only"
        ))
    
    return results


def fetch_rsid_info(
    rsids: List[str],
    rate_limit_s: float = 0.34,
    timeout_s: int = 15
) -> Dict[str, RsidInfo]:
    """
    Fetch information about rsIDs from Ensembl Variation API.
    
    Returns a dictionary mapping rsIDs to RsidInfo objects.
    """
    rsid_info_map: Dict[str, RsidInfo] = {}
    
    # Filter to only process rsIDs (start with "rs")
    valid_rsids = [r for r in rsids if r.startswith("rs")]
    
    for rsid in valid_rsids:
        url = f"https://rest.ensembl.org/variation/human/{rsid}"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        
        description = ""
        clinical_sig = ""
        phenotypes_str = ""
        minor_allele = ""
        maf = ""
        
        try:
            r = requests.get(url, headers=headers, timeout=timeout_s)
            if r.status_code == 200:
                data = r.json()
                
                # Extract clinical significance
                clinical_sigs = data.get("clinical_significance", [])
                if clinical_sigs:
                    clinical_sig = ", ".join(clinical_sigs)
                
                # Extract MAF info
                minor_allele = data.get("minor_allele", "") or ""
                maf_val = data.get("MAF")
                if maf_val is not None:
                    maf = f"{maf_val:.4f}"
                
                # Extract mappings for allele info
                mappings = data.get("mappings", [])
                if mappings:
                    m = mappings[0]
                    allele_str = m.get("allele_string", "")
                    location = m.get("location", "")
                    if allele_str and location:
                        description = f"{location} ({allele_str})"
                
                # Fetch phenotypes separately if available
                # Using phenotype endpoint
                pheno_url = f"https://rest.ensembl.org/variation/human/{rsid}?pops=1;phenotypes=1"
                try:
                    r2 = requests.get(pheno_url, headers=headers, timeout=timeout_s)
                    if r2.status_code == 200:
                        data2 = r2.json()
                        phenos = data2.get("phenotypes", [])
                        if phenos:
                            pheno_names = []
                            for p in phenos[:5]:  # Limit to first 5
                                trait = p.get("trait", "") or p.get("description", "")
                                if trait:
                                    pheno_names.append(trait)
                            if pheno_names:
                                phenotypes_str = "; ".join(pheno_names)
                except Exception:
                    pass
                
                rsid_info_map[rsid] = RsidInfo(
                    rsid=rsid,
                    description=description or rsid,
                    clinical_significance=clinical_sig or "Not reported",
                    phenotypes=phenotypes_str or "No phenotypes reported",
                    minor_allele=minor_allele or "N/A",
                    maf=maf or "N/A",
                    source="Ensembl Variation"
                )
                
                print(f"[API] rsID {rsid}: clinical={clinical_sig or 'N/A'}, phenotypes={'yes' if phenotypes_str else 'no'}")
        except Exception as e:
            rsid_info_map[rsid] = RsidInfo(
                rsid=rsid,
                description=rsid,
                clinical_significance="Fetch failed",
                phenotypes="Fetch failed",
                minor_allele="N/A",
                maf="N/A",
                source="Error"
            )
        
        time.sleep(rate_limit_s)
    
    return rsid_info_map


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
# Filtro de variantes centrais
# ----------------------------

def filter_central_variants(
    occs: List[VariantOcc], 
    central_window: Optional[int]
) -> Tuple[List[VariantOcc], List[VariantOcc]]:
    """
    Filtra ocorrências por distância ao centro.
    
    Args:
        occs: Lista de ocorrências de variantes
        central_window: Tamanho da janela central (em bp). 
                        Variantes dentro de ±(central_window/2) do centro são mantidas.
                        Se None, mantém todas.
    
    Returns:
        Tupla (occs_filtradas, occs_todas) onde:
        - occs_filtradas: apenas variantes dentro da janela central
        - occs_todas: todas as variantes (para comparação)
    """
    if central_window is None:
        return occs, occs
    
    half_window = central_window // 2
    filtered = [o for o in occs if o.distance_to_center <= half_window]
    
    return filtered, occs


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
    occs: List[VariantOcc],
    gene_info: Dict[str, GeneInfo],
    high_impact_info: Optional[List[HighImpactVariantInfo]] = None
) -> None:
    """
    Generate a textual report with:
      - overview
      - highlights by gene
      - highlights by individual
      - gene × individual matrix
    """
    n_occ = len(occs)
    uniq_vars = len({o.key for o in occs})

    genes = sorted(gene_sum.keys())
    indivs = sorted(indiv_sum.keys())

    # top genes by HIGH+MODERATE
    def score_gene(g: str) -> Tuple[int, int, int]:
        d = gene_sum[g]
        return (d.get("impact_HIGH", 0) + d.get("impact_MODERATE", 0),
                d.get("impact_HIGH", 0),
                d.get("unique_variants_count", 0))

    top_genes = sorted(genes, key=score_gene, reverse=True)[:15]

    # top individuals by HIGH+MODERATE
    def score_ind(i: str) -> Tuple[int, int, int]:
        d = indiv_sum[i]
        return (d.get("impact_HIGH", 0) + d.get("impact_MODERATE", 0),
                d.get("impact_HIGH", 0),
                d.get("unique_variants_count", 0))

    top_indivs = sorted(indivs, key=score_ind, reverse=True)[:15]

    with report_path.open("w", encoding="utf-8") as f:
        f.write("# Variant Annotation Report (by gene / by individual)\n\n")
        f.write("## Overview\n")
        f.write(f"- Occurrences (individual × window): **{n_occ}**\n")
        f.write(f"- Unique variants (chrom, pos, ref, alt): **{uniq_vars}**\n")
        f.write(f"- Genes (from input): **{len(genes)}**\n")
        f.write(f"- Individuals/hap (SAMPLE_HAP): **{len(indivs)}**\n\n")

        f.write("## Highlights by Gene (ordered by HIGH+MODERATE)\n\n")
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

        # Gene descriptions from API
        if gene_info:
            f.write("### Gene Descriptions (from APIs)\n\n")
            for g in genes:
                info = gene_info.get(g)
                if info and info.source != "No data found":
                    strand_info = f", strand={info.strand}" if info.strand and info.strand != "unknown" else ""
                    f.write(f"- **{g}** ({info.biotype}{strand_info}): {info.description}")
                    if info.function and info.function != "Function not available from APIs.":
                        truncated_func = info.function[:150] + "..." if len(info.function) > 150 else info.function
                        f.write(f" Function: {truncated_func}")
                    f.write(f" [Source: {info.source}]\n")
            f.write("\n")

        f.write("## Highlights by Individual (SAMPLE_HAP) (ordered by HIGH+MODERATE)\n\n")
        f.write("| Individual | Unique variants | Occurrences | HIGH | MODERATE | LOW | MODIFIER | Top consequences | rsIDs |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---|---|\n")
        for i_id in top_indivs:
            d = indiv_sum[i_id]
            f.write(
                f"| {i_id} | {d['unique_variants_count']} | {d['total_occurrences']} | "
                f"{d.get('impact_HIGH',0)} | {d.get('impact_MODERATE',0)} | {d.get('impact_LOW',0)} | {d.get('impact_MODIFIER',0)} | "
                f"{d.get('top_msc_str','')} | {d.get('rsids_str','')} |\n"
            )
        f.write("\n")

        f.write("## Individual × Gene Matrix (summary)\n\n")
        f.write("Below are the pairs (individual, gene) with highest HIGH+MODERATE count.\n\n")
        pairs = list(indiv_gene_sum.keys())

        def score_pair(p: Tuple[str, str]) -> Tuple[int, int, int]:
            d = indiv_gene_sum[p]
            return (d.get("impact_HIGH", 0) + d.get("impact_MODERATE", 0),
                    d.get("impact_HIGH", 0),
                    d.get("unique_variants_count", 0))

        top_pairs = sorted(pairs, key=score_pair, reverse=True)[:30]
        f.write("| Individual | Gene | Unique variants | Occurrences | HIGH | MODERATE | Top consequences | rsIDs |\n")
        f.write("|---|---|---:|---:|---:|---:|---|---|\n")
        for (i_id, g) in top_pairs:
            d = indiv_gene_sum[(i_id, g)]
            f.write(
                f"| {i_id} | {g} | {d['unique_variants_count']} | {d['total_occurrences']} | "
                f"{d.get('impact_HIGH',0)} | {d.get('impact_MODERATE',0)} | {d.get('top_msc_str','')} | {d.get('rsids_str','')} |\n"
            )
        f.write("\n")

        # Interpretation notes
        # HIGH Impact Variants - Detailed Analysis
        if high_impact_info and len(high_impact_info) > 0:
            f.write("## HIGH Impact Variants - Detailed Analysis\n\n")
            f.write(f"Found **{len(high_impact_info)} HIGH impact variant(s)** that may significantly affect protein function:\n\n")
            
            for i, hv in enumerate(high_impact_info, 1):
                f.write(f"### {i}. {hv.chrom}:{hv.pos} {hv.ref}>{hv.alt} ({hv.gene})\n\n")
                f.write(f"| Property | Value |\n")
                f.write(f"|----------|-------|\n")
                f.write(f"| **Gene** | {hv.gene} |\n")
                f.write(f"| **Consequence** | {hv.consequence} |\n")
                if hv.amino_acids:
                    f.write(f"| **Amino acid change** | {hv.amino_acids} |\n")
                if hv.codons:
                    f.write(f"| **Codon change** | {hv.codons} |\n")
                f.write(f"| **Distance to DeepLIFT center** | {hv.distance_to_center}bp |\n")
                if hv.rsid:
                    f.write(f"| **rsID** | [{hv.rsid}]({hv.ensembl_link}) |\n")
                if hv.cosmic_id:
                    f.write(f"| **COSMIC ID** | {hv.cosmic_id} |\n")
                f.write(f"| **Clinical significance** | {hv.clinical_significance} |\n")
                f.write(f"| **Associated phenotypes** | {hv.phenotypes} |\n")
                if hv.maf != "N/A":
                    f.write(f"| **Minor allele frequency** | {hv.maf} |\n")
                if hv.transcript_id:
                    f.write(f"| **Transcript** | {hv.transcript_id} |\n")
                f.write(f"| **Data source** | {hv.source} |\n")
                f.write("\n")
                
                # Add interpretation for specific consequence types
                if hv.consequence == "stop_gained":
                    f.write(f"> ⚠️ **Stop gained**: This variant introduces a premature stop codon, ")
                    f.write(f"likely resulting in a truncated, non-functional protein or nonsense-mediated decay (NMD).\n\n")
                elif hv.consequence == "splice_donor_variant":
                    f.write(f"> ⚠️ **Splice donor disruption**: This variant affects the splice donor site, ")
                    f.write(f"potentially causing exon skipping or intron retention, leading to abnormal protein.\n\n")
                elif hv.consequence == "splice_acceptor_variant":
                    f.write(f"> ⚠️ **Splice acceptor disruption**: This variant affects the splice acceptor site, ")
                    f.write(f"potentially causing aberrant splicing and abnormal protein production.\n\n")
                elif hv.consequence == "frameshift_variant":
                    f.write(f"> ⚠️ **Frameshift**: This insertion/deletion shifts the reading frame, ")
                    f.write(f"likely producing a completely different and non-functional protein.\n\n")
            
            f.write("---\n\n")
        
        f.write("## Interpretation Notes\n\n")
        f.write("- **VEP** indicates **functional consequences** (e.g., missense, splice, intronic) and often returns **rsIDs** when the variant matches a cataloged one.\n")
        f.write("- This helps prioritize candidates, but does **not prove** increased/decreased expression without additional data (eQTL/expression/proteomics).\n")
        f.write("- This report **highlights** genes and individuals with (i) higher variant burden, (ii) more severe consequences, and (iii) presence of known rsIDs.\n\n")

        f.write("## Practical Notes\n")
        f.write(
            '1) To infer "increased protein production", prioritize variants in **coding regions** (missense/nonsense/frameshift) and **splicing**; '
            "for expression effects, cross-reference with eQTL/regulatory tracks.\n"
            "2) If VEP returns rsIDs (colocated_variants), you can cross-reference with literature/GTEx/ClinVar/Ensembl Variation.\n"
            "3) For large batches, local VEP (script + cache) tends to be more stable than REST API.\n"
        )


def write_phenotype_validation_report(
    report_path: Path,
    occs: List[VariantOcc],
    occs_all: List[VariantOcc],
    ann: Dict[VariantKey, VariantAnnotation],
    central_window: Optional[int],
    window_size: int,
    gene_sum: Dict[str, dict],
    gene_info: Dict[str, GeneInfo],
    rsid_info: Dict[str, RsidInfo]
) -> None:
    """
    Generate a generic phenotype validation report.
    
    This report focuses on:
    1. Variants near the center identified by DeepLIFT
    2. Genes detected in the data (with information from APIs)
    3. rsIDs found in the variants (with information from Ensembl Variation)
    4. Statistical summary and next steps
    """
    n_filtered = len(occs)
    n_all = len(occs_all)
    unique_filtered = len({o.key for o in occs})
    unique_all = len({o.key for o in occs_all})
    
    # Get detected genes from data
    detected_genes = sorted({o.gene for o in occs})
    
    # Collect all rsIDs from annotations
    all_rsids_in_data: Dict[str, List[Tuple[VariantKey, VariantAnnotation, int]]] = {}
    for o in occs:
        ann_entry = ann.get(o.key)
        if ann_entry and ann_entry.rsids:
            for rsid in ann_entry.rsids.split(","):
                rsid = rsid.strip()
                if rsid and rsid.startswith("rs"):
                    if rsid not in all_rsids_in_data:
                        all_rsids_in_data[rsid] = []
                    all_rsids_in_data[rsid].append((o.key, ann_entry, o.distance_to_center))
    
    # Count variants by impact
    total_high = 0
    total_moderate = 0
    total_missense = 0
    for o in occs:
        ann_entry = ann.get(o.key)
        if ann_entry:
            if ann_entry.worst_impact == "HIGH":
                total_high += 1
            elif ann_entry.worst_impact == "MODERATE":
                total_moderate += 1
            if "missense_variant" in ann_entry.consequence_terms:
                total_missense += 1
    
    # Calculate distance statistics
    if occs:
        distances = [o.distance_to_center for o in occs]
        avg_distance = sum(distances) / len(distances)
        min_distance = min(distances)
        max_distance = max(distances)
    else:
        avg_distance = min_distance = max_distance = 0
    
    with report_path.open("w", encoding="utf-8") as f:
        f.write("# Phenotype Validation Report\n\n")
        
        f.write("## Overview\n\n")
        f.write("This report analyzes variants identified by the DeepLIFT → VEP pipeline ")
        f.write("to understand which genetic regions contribute to the phenotype of interest.\n\n")
        
        f.write("---\n\n")
        
        f.write("## Statistical Summary\n\n")
        if central_window is not None:
            half_win = central_window // 2
            f.write("| Metric | Value |\n")
            f.write("|--------|-------|\n")
            f.write(f"| Central window | ±{half_win}bp from DeepLIFT center |\n")
            f.write(f"| Total variants ({window_size}bp) | {n_all} occurrences, {unique_all} unique |\n")
            f.write(f"| Filtered variants (central window) | {n_filtered} occurrences, {unique_filtered} unique |\n")
            if n_all > 0:
                f.write(f"| Reduction | {100*(1 - n_filtered/n_all):.1f}% |\n")
            f.write(f"| Mean distance to center | {avg_distance:.1f}bp |\n")
            f.write(f"| Distance range | {min_distance}-{max_distance}bp |\n")
        else:
            f.write("| Metric | Value |\n")
            f.write("|--------|-------|\n")
            f.write(f"| Window | {window_size}bp (no central filter) |\n")
            f.write(f"| Total variants | {n_all} occurrences, {unique_all} unique |\n")
            f.write(f"| Mean distance to center | {avg_distance:.1f}bp |\n")
        f.write("\n")
        
        f.write("---\n\n")
        
        f.write("## Detected Genes\n\n")
        f.write(f"The following **{len(detected_genes)} genes** were identified in regions with high DeepLIFT attribution:\n\n")
        
        for gene in detected_genes:
            gene_data = gene_sum.get(gene, {})
            info = gene_info.get(gene)
            
            f.write(f"### {gene}\n\n")
            
            if info and info.source != "No data found":
                f.write(f"- **Description**: {info.description}\n")
                f.write(f"- **Function**: {info.function}\n")
                f.write(f"- **Biotype**: {info.biotype}\n")
                f.write(f"- **Chromosome**: {info.chromosome}\n")
                f.write(f"- **Expression strand**: {info.strand} (protein coding direction)\n")
                f.write(f"- **Data source**: {info.source}\n")
            else:
                f.write(f"- **Description**: No API data available for this gene.\n")
            
            f.write(f"- **Unique variants**: {gene_data.get('unique_variants_count', 0)}\n")
            f.write(f"- **HIGH impact**: {gene_data.get('impact_HIGH', 0)}\n")
            f.write(f"- **MODERATE impact**: {gene_data.get('impact_MODERATE', 0)}\n")
            f.write(f"- **Top consequences**: {gene_data.get('top_msc_str', 'N/A')}\n\n")
        
        f.write("---\n\n")
        
        f.write("## Annotated rsIDs\n\n")
        
        # Filter to rsIDs that we have info for
        rsids_with_info = {rsid: data for rsid, data in all_rsids_in_data.items() 
                          if rsid in rsid_info}
        rsids_without_info = {rsid: data for rsid, data in all_rsids_in_data.items() 
                             if rsid not in rsid_info}
        
        if rsids_with_info:
            f.write(f"**{len(rsids_with_info)} rsIDs** found in the central window with Ensembl Variation annotations:\n\n")
            
            # Sort by number of occurrences
            sorted_rsids = sorted(rsids_with_info.items(), 
                                  key=lambda x: (-len(x[1]), x[0]))[:20]  # Top 20
            
            for rsid, occurrences in sorted_rsids:
                info = rsid_info[rsid]
                min_dist = min(d for _, _, d in occurrences)
                key, ann_entry, _ = occurrences[0]
                
                f.write(f"### {rsid}\n\n")
                f.write(f"- **Position**: {info.description}\n")
                f.write(f"- **Clinical significance**: {info.clinical_significance}\n")
                f.write(f"- **Associated phenotypes**: {info.phenotypes}\n")
                f.write(f"- **Minor allele**: {info.minor_allele} (MAF: {info.maf})\n")
                f.write(f"- **Occurrences in data**: {len(occurrences)}\n")
                f.write(f"- **Min distance to DeepLIFT center**: {min_dist}bp\n")
                f.write(f"- **Consequence**: {ann_entry.most_severe_consequence}\n\n")
            
            if len(rsids_with_info) > 20:
                f.write(f"*... and {len(rsids_with_info) - 20} more rsIDs*\n\n")
        else:
            f.write("No rsIDs with Ensembl Variation annotations were found in the central window.\n\n")
        
        if rsids_without_info:
            f.write(f"\n**{len(rsids_without_info)} additional rsIDs** found without detailed annotations:\n")
            f.write(", ".join(sorted(rsids_without_info.keys())[:30]))
            if len(rsids_without_info) > 30:
                f.write(f" ... and {len(rsids_without_info) - 30} more")
            f.write("\n\n")
        
        f.write("---\n\n")
        
        f.write("## Functional Impact Analysis\n\n")
        f.write("Variants with potential functional impact:\n\n")
        f.write("| Category | Count |\n")
        f.write("|----------|-------|\n")
        f.write(f"| HIGH impact (stop, frameshift, splice) | {total_high} |\n")
        f.write(f"| MODERATE impact (missense, in-frame) | {total_moderate} |\n")
        f.write(f"| Missense variants | {total_missense} |\n")
        f.write("\n")
        
        f.write("---\n\n")
        
        f.write("## Summary\n\n")
        
        # Evaluate evidence
        evidence_points = []
        
        if detected_genes:
            evidence_points.append(f"- {len(detected_genes)} gene(s) identified by DeepLIFT")
        
        if rsids_with_info:
            # Check for clinically significant variants
            clinically_significant = [r for r, info in rsid_info.items() 
                                      if info.clinical_significance not in ("Not reported", "Fetch failed", "")]
            if clinically_significant:
                evidence_points.append(f"- {len(clinically_significant)} rsID(s) with clinical significance annotations")
            
            # Check for phenotype associations
            with_phenotypes = [r for r, info in rsid_info.items() 
                              if info.phenotypes not in ("No phenotypes reported", "Fetch failed", "")]
            if with_phenotypes:
                evidence_points.append(f"- {len(with_phenotypes)} rsID(s) with associated phenotypes")
        
        if total_high > 0:
            evidence_points.append(f"- {total_high} HIGH impact variant(s) detected")
        
        if total_moderate > 0:
            evidence_points.append(f"- {total_moderate} MODERATE impact variant(s) detected")
        
        if central_window and n_filtered > 0:
            evidence_points.append(f"- {n_filtered} variant(s) in the central region (DeepLIFT focus)")
        
        if evidence_points:
            f.write("### Key Findings\n\n")
            for point in evidence_points:
                f.write(f"{point}\n")
            f.write("\n")
        else:
            f.write("No significant findings in the analyzed region.\n\n")
        
        f.write("### Recommended Next Steps\n\n")
        f.write("1. **Cross-reference with GWAS**: Compare detected variants with published GWAS for related phenotypes\n")
        f.write("2. **Check eQTL databases**: Verify if variants affect gene expression (GTEx, eQTLGen)\n")
        f.write("3. **Functional validation**: Consider experimental validation for HIGH/MODERATE impact variants\n")
        f.write("4. **Population analysis**: Check allele frequencies across populations (gnomAD, 1000 Genomes)\n")
        f.write("5. **Literature review**: Search for gene-phenotype associations in PubMed/OMIM\n\n")
        
        f.write("---\n\n")
        f.write("*Report generated automatically by annotate_deeplift_windows.py*\n")


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
    ap.add_argument(
        "--haplotype",
        type=str,
        choices=["H1", "H2", "both"],
        default="H1",
        help="Qual haplótipo considerar: H1 (default), H2, ou both (ambos)."
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
    validation_report_path = outdir / "phenotype_validation_report.md"

    records = parse_input_file(input_path)

    # Filtrar por haplótipo
    total_before_hap_filter = len(records)
    if args.haplotype == "H1":
        records = [r for r in records if r.hap.upper() == "H1"]
    elif args.haplotype == "H2":
        records = [r for r in records if r.hap.upper() == "H2"]
    # else: "both" - mantém todos
    
    total_after_hap_filter = len(records)
    if args.haplotype != "both":
        print(f"[INFO] Haplótipo: {args.haplotype} (filtrados {total_after_hap_filter} de {total_before_hap_filter} records)")
    else:
        print(f"[INFO] Haplótipo: ambos (H1 + H2), {total_after_hap_filter} records")
    
    if not records:
        print(f"[ERRO] Nenhum record encontrado para haplótipo {args.haplotype}. Encerrando.")
        return 1

    # Detectar tamanho da janela a partir dos dados
    window_size = len(records[0].seq)
    seq_sizes = {len(r.seq) for r in records}
    if len(seq_sizes) > 1:
        print(f"[WARN] Sequências com tamanhos diferentes detectadas: {sorted(seq_sizes)}")
        print(f"[WARN] Usando tamanho da primeira sequência: {window_size}bp")
    print(f"[INFO] Tamanho da janela detectado: {window_size}bp")

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
        start0, end1 = window_centered(rec.center_1based, window_size)
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
            # Calcular distância ao centro da janela
            distance = abs(k.pos_1based - rec.center_1based)
            occs.append(VariantOcc(
                key=k, 
                sample=rec.sample, 
                hap=rec.hap, 
                gene=rec.gene, 
                header=rec.header,
                center_1based=rec.center_1based,
                distance_to_center=distance
            ))
            unique_keys_set[k] = None

    if not occs:
        print("[INFO] Nenhuma variante detectada. Encerrando.")
        return 0

    # Aplicar filtro central se especificado
    occs_all = occs  # Manter referência para estatísticas
    occs, occs_all = filter_central_variants(occs, args.central_window)
    
    # Recalcular unique_keys apenas para variantes filtradas
    unique_keys_set_filtered: Dict[VariantKey, None] = {}
    for o in occs:
        unique_keys_set_filtered[o.key] = None
    unique_keys = list(unique_keys_set_filtered.keys())
    
    # Estatísticas
    total_all = len(occs_all)
    total_filtered = len(occs)
    unique_all = len({o.key for o in occs_all})
    unique_filtered = len(unique_keys)
    
    if args.central_window is not None:
        half_win = args.central_window // 2
        print(f"[INFO] Filtro central: ±{half_win}bp do centro (janela de {args.central_window}bp)")
        print(f"[INFO] Variantes originais: {total_all} ocorrências, {unique_all} únicas")
        print(f"[INFO] Variantes filtradas: {total_filtered} ocorrências, {unique_filtered} únicas")
        print(f"[INFO] Redução: {100*(1 - total_filtered/total_all):.1f}% das ocorrências")
    else:
        print(f"[INFO] Variant occurrences: {total_filtered} | unique variants: {unique_filtered}")

    # Salvar ocorrências (com distância ao centro)
    occ_header = ["chrom", "pos_1based", "ref", "alt", "distance_to_center", "center_1based", "sample", "hap", "gene", "header"]
    occ_rows = [
        [o.key.chrom, str(o.key.pos_1based), o.key.ref, o.key.alt, 
         str(o.distance_to_center), str(o.center_1based),
         o.sample, o.hap, o.gene, o.header] 
        for o in occs
    ]
    write_tsv(occs_path, occ_header, occ_rows)
    print(f"[SAVED] variant occurrences TSV: {occs_path}")

    # Calcular distância mínima para cada variante única (a mesma variante pode aparecer em múltiplos indivíduos)
    variant_min_distance: Dict[VariantKey, int] = {}
    for o in occs:
        if o.key not in variant_min_distance:
            variant_min_distance[o.key] = o.distance_to_center
        else:
            variant_min_distance[o.key] = min(variant_min_distance[o.key], o.distance_to_center)
    
    # Salvar variantes únicas (com distância mínima ao centro)
    uv_header = ["chrom", "pos_1based", "ref", "alt", "min_distance_to_center"]
    uv_rows = [
        [k.chrom, str(k.pos_1based), k.ref, k.alt, str(variant_min_distance.get(k, 0))] 
        for k in sorted(unique_keys, key=lambda x: (x.chrom, x.pos_1based, x.ref, x.alt))
    ]
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

        # Salvar anotações (TSV) - inclui distância mínima ao centro
        va_header = [
            "chrom", "pos_1based", "ref", "alt", "min_distance_to_center",
            "most_severe_consequence", "worst_impact",
            "gene_symbols", "gene_ids", "transcript_ids", "protein_ids",
            "consequence_terms", "amino_acids", "codons", "rsids"
        ]
        va_rows: List[List[str]] = []
        for k, a in sorted(ann_map.items(), key=lambda x: (x[0].chrom, x[0].pos_1based, x[0].ref, x[0].alt)):
            va_rows.append([
                k.chrom, str(k.pos_1based), k.ref, k.alt, 
                str(variant_min_distance.get(k, 0)),
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

    # --- Fetch gene information from APIs ---
    detected_genes = sorted({o.gene for o in occs})
    print(f"[INFO] Fetching information for {len(detected_genes)} genes from APIs...")
    gene_info = fetch_gene_info(detected_genes)
    print(f"[INFO] Gene information fetched for {len([g for g in gene_info.values() if g.source != 'No data found'])} genes")

    # --- Fetch rsID information from APIs ---
    # Collect all rsIDs from VEP annotations
    all_rsids: List[str] = []
    for a in ann_map.values():
        if a.rsids:
            for rsid in a.rsids.split(","):
                rsid = rsid.strip()
                if rsid and rsid.startswith("rs"):
                    all_rsids.append(rsid)
    unique_rsids = sorted(set(all_rsids))
    
    # Limit to first 50 rsIDs to avoid too many API calls
    rsids_to_fetch = unique_rsids[:50]
    rsid_info: Dict[str, RsidInfo] = {}
    if rsids_to_fetch:
        print(f"[INFO] Fetching information for {len(rsids_to_fetch)} rsIDs from Ensembl Variation...")
        rsid_info = fetch_rsid_info(rsids_to_fetch)
        print(f"[INFO] rsID information fetched for {len(rsid_info)} variants")
    else:
        print("[INFO] No rsIDs to fetch (VEP may be disabled or no colocated variants)")

    # Fetch detailed information for HIGH impact variants
    high_impact_info: List[HighImpactVariantInfo] = []
    high_variants = [
        (k, a, variant_min_distance.get(k, 0)) 
        for k, a in ann_map.items() 
        if a.worst_impact == "HIGH"
    ]
    if high_variants:
        print(f"[INFO] Fetching detailed information for {len(high_variants)} HIGH impact variant(s)...")
        high_impact_info = fetch_high_impact_details(high_variants)
        print(f"[INFO] HIGH impact details fetched for {len(high_impact_info)} variant(s)")

    # relatório textual
    write_report_markdown(report_path, gene_sum, indiv_sum, indiv_gene_sum, ann_map, occs, gene_info, high_impact_info)
    print(f"[SAVED] report (Markdown): {report_path}")

    # Relatório de validação de fenótipo
    write_phenotype_validation_report(
        validation_report_path, 
        occs, 
        occs_all, 
        ann_map, 
        args.central_window,
        window_size,
        gene_sum,
        gene_info,
        rsid_info
    )
    print(f"[SAVED] phenotype validation report (Markdown): {validation_report_path}")

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
    print(f"[OUT] {validation_report_path}")
    print("==========================\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

