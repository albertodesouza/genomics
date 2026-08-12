#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
validate_pigmentation_hypothesis.py

Valida a hipótese de que o pipeline DeepLIFT → VEP identifica variantes
genéticas associadas à pigmentação em africanos.

Três níveis de validação:
  1. Genes: Verifica se os genes detectados correspondem a genes de pigmentação conhecidos
  2. Variantes: Verifica se os rsIDs encontrados incluem SNPs conhecidos de GWAS
  3. Mecanismo: Analisa se as variantes de alto impacto fazem sentido biológico

Uso:
  python3 validate_pigmentation_hypothesis.py <input_dir> [--output <output_file>]

Exemplo:
  python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 --output validation_report.md
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Optional: requests for gnomAD API
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False


# ============================================================================
# GENE DATABASE - Genes de pigmentação conhecidos com referências GWAS
# ============================================================================

@dataclass
class PigmentationGene:
    """Informação sobre um gene de pigmentação conhecido."""
    name: str
    category: str  # "primary", "regulatory", "uv_response", "other"
    function: str
    omim: str
    gwas_references: List[str]
    population_effect: str
    known_rsids: List[str] = field(default_factory=list)


# Banco de dados de genes de pigmentação conhecidos
PIGMENTATION_GENES_DB: Dict[str, PigmentationGene] = {
    "SLC24A5": PigmentationGene(
        name="SLC24A5",
        category="primary",
        function="Calcium antiporter in melanosomes, critical for melanin synthesis",
        omim="609802",
        gwas_references=["Lamason et al. 2005 Science", "Crawford et al. 2017 Science"],
        population_effect="Light skin in Europeans (rs1426654 A allele)",
        known_rsids=["rs1426654"]
    ),
    "SLC45A2": PigmentationGene(
        name="SLC45A2",
        category="primary",
        function="Melanosomal membrane transport protein, melanin synthesis",
        omim="606202",
        gwas_references=["Graf et al. 2005", "Stokowski et al. 2007"],
        population_effect="Light skin/hair in Europeans (rs16891982)",
        known_rsids=["rs16891982", "rs26722"]
    ),
    "OCA2": PigmentationGene(
        name="OCA2",
        category="primary",
        function="Melanosomal transmembrane protein, modulates pH for tyrosinase activity",
        omim="611409",
        gwas_references=["Sturm et al. 2008", "Eiberg et al. 2008"],
        population_effect="Eye color, skin pigmentation, oculocutaneous albinism type 2",
        known_rsids=["rs1800407", "rs12913832"]  # rs12913832 is actually in HERC2 but affects OCA2
    ),
    "HERC2": PigmentationGene(
        name="HERC2",
        category="regulatory",
        function="E3 ubiquitin ligase, regulates OCA2 expression",
        omim="605837",
        gwas_references=["Sturm et al. 2008", "Eiberg et al. 2008"],
        population_effect="Blue/brown eye color (rs12913832)",
        known_rsids=["rs12913832", "rs916977"]
    ),
    "TYR": PigmentationGene(
        name="TYR",
        category="primary",
        function="Tyrosinase, rate-limiting enzyme in melanin biosynthesis",
        omim="606933",
        gwas_references=["Sulem et al. 2007", "Han et al. 2008"],
        population_effect="Oculocutaneous albinism type 1, skin/hair/eye color",
        known_rsids=["rs1042602", "rs1126809"]
    ),
    "TYRP1": PigmentationGene(
        name="TYRP1",
        category="primary",
        function="Tyrosinase-related protein 1, melanin biosynthesis stabilization",
        omim="115501",
        gwas_references=["Sulem et al. 2007"],
        population_effect="Hair/skin color variation",
        known_rsids=["rs1408799", "rs683"]
    ),
    "MC1R": PigmentationGene(
        name="MC1R",
        category="primary",
        function="Melanocortin 1 receptor, switches eumelanin/pheomelanin production",
        omim="155555",
        gwas_references=["Valverde et al. 1995", "Box et al. 1997"],
        population_effect="Red hair, fair skin, freckling",
        known_rsids=["rs1805007", "rs1805008", "rs1805009"]
    ),
    "ASIP": PigmentationGene(
        name="ASIP",
        category="regulatory",
        function="Agouti signaling protein, MC1R antagonist",
        omim="600201",
        gwas_references=["Kanetsky et al. 2002"],
        population_effect="Skin darkness variation",
        known_rsids=["rs4911414", "rs1015362"]
    ),
    "DDB1": PigmentationGene(
        name="DDB1",
        category="uv_response",
        function="DNA damage binding protein 1, UV damage repair",
        omim="600045",
        gwas_references=["Crawford et al. 2017 Science"],
        population_effect="UV damage response, identified in African pigmentation GWAS",
        known_rsids=["rs11230664", "rs7120594"]
    ),
    "MFSD12": PigmentationGene(
        name="MFSD12",
        category="primary",
        function="Major facilitator superfamily domain containing 12, lysosomal function",
        omim="618170",
        gwas_references=["Crawford et al. 2017 Science"],
        population_effect="Skin pigmentation in African populations",
        known_rsids=["rs10424065", "rs56203814"]
    ),
    "EDAR": PigmentationGene(
        name="EDAR",
        category="other",
        function="Ectodysplasin A receptor, hair morphology and sweat gland density",
        omim="604095",
        gwas_references=["Fujimoto et al. 2008", "Mou et al. 2008"],
        population_effect="Hair thickness, shovel-shaped incisors (East Asian)",
        known_rsids=["rs3827760"]
    ),
    "KITLG": PigmentationGene(
        name="KITLG",
        category="other",
        function="KIT ligand, melanocyte development and survival",
        omim="184745",
        gwas_references=["Miller et al. 2007", "Sulem et al. 2007"],
        population_effect="Hair color variation",
        known_rsids=["rs12821256"]
    ),
    "IRF4": PigmentationGene(
        name="IRF4",
        category="regulatory",
        function="Interferon regulatory factor 4, regulates TYR expression",
        omim="601900",
        gwas_references=["Han et al. 2008", "Praetorius et al. 2013"],
        population_effect="Skin sensitivity to sun, freckling, hair color",
        known_rsids=["rs12203592"]
    ),
    "BNC2": PigmentationGene(
        name="BNC2",
        category="regulatory",
        function="Basonuclin 2, transcription factor in melanocytes",
        omim="608669",
        gwas_references=["Jacobs et al. 2013"],
        population_effect="Skin saturation, freckling",
        known_rsids=["rs2153271", "rs10756819"]
    ),
}

# SNPs conhecidos de pigmentação com seus efeitos
KNOWN_PIGMENTATION_RSIDS: Dict[str, Dict] = {
    "rs1426654": {
        "gene": "SLC24A5",
        "effect": "Light skin in Europeans (A allele: ~100% European, ~0% African)",
        "consequence": "missense_variant (Ala111Thr)",
        "references": ["Lamason et al. 2005 Science"]
    },
    "rs16891982": {
        "gene": "SLC45A2",
        "effect": "Light skin/hair in Europeans (G allele ancestral, C derived)",
        "consequence": "missense_variant (Leu374Phe)",
        "references": ["Graf et al. 2005", "Stokowski et al. 2007"]
    },
    "rs12913832": {
        "gene": "HERC2 (regulates OCA2)",
        "effect": "Blue/brown eye color (A: brown, G: blue)",
        "consequence": "intron_variant (regulatory)",
        "references": ["Sturm et al. 2008", "Eiberg et al. 2008"]
    },
    "rs1800407": {
        "gene": "OCA2",
        "effect": "Eye color modifier (Arg419Gln)",
        "consequence": "missense_variant",
        "references": ["Rebbeck et al. 2002"]
    },
    "rs1042602": {
        "gene": "TYR",
        "effect": "Skin/hair color (Ser192Tyr)",
        "consequence": "missense_variant",
        "references": ["Shriver et al. 2003"]
    },
    "rs1805007": {
        "gene": "MC1R",
        "effect": "Red hair, fair skin (Arg151Cys - R allele)",
        "consequence": "missense_variant",
        "references": ["Valverde et al. 1995"]
    },
    "rs1805008": {
        "gene": "MC1R",
        "effect": "Red hair, fair skin (Arg160Trp - r allele)",
        "consequence": "missense_variant",
        "references": ["Valverde et al. 1995"]
    },
    "rs11230664": {
        "gene": "DDB1",
        "effect": "Associated with darker skin in African populations",
        "consequence": "intron_variant",
        "references": ["Crawford et al. 2017 Science"]
    },
    "rs10424065": {
        "gene": "MFSD12",
        "effect": "African-specific skin pigmentation variant",
        "consequence": "regulatory_region_variant",
        "references": ["Crawford et al. 2017 Science"]
    },
    "rs3827760": {
        "gene": "EDAR",
        "effect": "Hair thickness, sweat glands (Val370Ala, East Asian)",
        "consequence": "missense_variant",
        "references": ["Fujimoto et al. 2008"]
    },
    "rs12203592": {
        "gene": "IRF4",
        "effect": "Freckling, sun sensitivity, hair color",
        "consequence": "intron_variant (regulatory)",
        "references": ["Han et al. 2008"]
    },
}

# Crawford et al. 2017 (Science) - landmark African pigmentation GWAS genes
CRAWFORD_2017_GENES = {"SLC24A5", "MFSD12", "DDB1", "HERC2", "OCA2"}


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class DetectedGene:
    """Gene detectado pelo pipeline."""
    name: str
    unique_variants: int
    high_impact: int
    moderate_impact: int
    low_impact: int
    modifier_impact: int
    top_consequences: str
    rsids: List[str]


@dataclass
class VariantInfo:
    """Informação de uma variante do VEP."""
    chrom: str
    pos: int
    ref: str
    alt: str
    distance_to_center: int
    consequence: str
    impact: str
    gene: str
    gene_id: str
    transcript_ids: str
    protein_ids: str
    all_consequences: str
    amino_acids: str
    codons: str
    rsids: str


@dataclass
class ValidationResult:
    """Resultado da validação."""
    # Level 1: Genes
    detected_genes: List[str]
    matched_pigmentation_genes: List[str]
    unmatched_genes: List[str]
    gene_match_percentage: float
    primary_genes_count: int
    regulatory_genes_count: int
    uv_response_genes_count: int
    crawford_overlap_count: int
    crawford_overlap_percentage: float
    
    # Level 2: Variants
    all_rsids_found: List[str]
    known_pigmentation_rsids_found: List[str]
    total_high_impact: int
    total_moderate_impact: int
    total_missense: int
    
    # Level 3: Mechanism
    mechanism_summary: Dict[str, str]
    
    # Population frequencies (optional)
    population_frequencies: Optional[Dict[str, Dict]] = None
    
    # Overall score
    validation_status: str = ""
    validation_score: float = 0.0


# ============================================================================
# PARSING FUNCTIONS
# ============================================================================

def parse_summary_by_gene(filepath: Path) -> Dict[str, DetectedGene]:
    """Parse summary_by_gene.tsv file."""
    genes = {}
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rsids = [r.strip() for r in row.get("rsids", "").split(",") if r.strip()]
            genes[row["gene"]] = DetectedGene(
                name=row["gene"],
                unique_variants=int(row.get("unique_variants", 0)),
                high_impact=int(row.get("HIGH", 0)),
                moderate_impact=int(row.get("MODERATE", 0)),
                low_impact=int(row.get("LOW", 0)),
                modifier_impact=int(row.get("MODIFIER", 0)),
                top_consequences=row.get("top_consequences", ""),
                rsids=rsids
            )
    return genes


def parse_vep_annotated_variants(filepath: Path) -> List[VariantInfo]:
    """Parse vep_annotated_variants.tsv file."""
    variants = []
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variants.append(VariantInfo(
                chrom=row.get("chrom", ""),
                pos=int(row.get("pos_1based", row.get("pos", 0))),
                ref=row.get("ref", ""),
                alt=row.get("alt", ""),
                distance_to_center=int(row.get("min_distance_to_center", row.get("distance_to_center", 0))),
                consequence=row.get("most_severe_consequence", ""),
                impact=row.get("worst_impact", ""),
                gene=row.get("gene_symbols", ""),
                gene_id=row.get("gene_id", ""),
                transcript_ids=row.get("transcript_ids", ""),
                protein_ids=row.get("protein_ids", ""),
                all_consequences=row.get("consequence_terms", ""),
                amino_acids=row.get("amino_acids", ""),
                codons=row.get("codons", ""),
                rsids=row.get("rsids", "")
            ))
    return variants


# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

def validate_genes(detected_genes: Dict[str, DetectedGene]) -> Tuple[List[str], List[str], int, int, int, int, float]:
    """
    Level 1: Validate detected genes against known pigmentation genes.
    
    Returns: (matched, unmatched, primary_count, regulatory_count, uv_count, crawford_overlap, crawford_percentage)
    """
    matched = []
    unmatched = []
    primary_count = 0
    regulatory_count = 0
    uv_count = 0
    
    for gene_name in detected_genes:
        if gene_name in PIGMENTATION_GENES_DB:
            matched.append(gene_name)
            gene_info = PIGMENTATION_GENES_DB[gene_name]
            if gene_info.category == "primary":
                primary_count += 1
            elif gene_info.category == "regulatory":
                regulatory_count += 1
            elif gene_info.category == "uv_response":
                uv_count += 1
        else:
            unmatched.append(gene_name)
    
    # Crawford et al. 2017 overlap
    detected_set = set(detected_genes.keys())
    crawford_overlap = len(detected_set & CRAWFORD_2017_GENES)
    crawford_percentage = crawford_overlap / len(CRAWFORD_2017_GENES) * 100 if CRAWFORD_2017_GENES else 0
    
    return matched, unmatched, primary_count, regulatory_count, uv_count, crawford_overlap, crawford_percentage


def validate_variants(detected_genes: Dict[str, DetectedGene], variants: List[VariantInfo]) -> Tuple[List[str], List[str], int, int, int]:
    """
    Level 2: Validate detected variants/rsIDs against known pigmentation SNPs.
    
    Returns: (all_rsids, known_pigmentation_rsids, high_count, moderate_count, missense_count)
    """
    # Collect all rsIDs
    all_rsids: Set[str] = set()
    for gene in detected_genes.values():
        for rsid in gene.rsids:
            if rsid.startswith("rs"):
                all_rsids.add(rsid)
    
    # Check which are known pigmentation rsIDs
    known_found = [rsid for rsid in all_rsids if rsid in KNOWN_PIGMENTATION_RSIDS]
    
    # Count impact levels and missense
    high_count = sum(g.high_impact for g in detected_genes.values())
    moderate_count = sum(g.moderate_impact for g in detected_genes.values())
    
    missense_count = 0
    for v in variants:
        if "missense_variant" in v.consequence or "missense_variant" in v.all_consequences:
            missense_count += 1
    
    return list(all_rsids), known_found, high_count, moderate_count, missense_count


def validate_mechanism(detected_genes: Dict[str, DetectedGene], variants: List[VariantInfo]) -> Dict[str, str]:
    """
    Level 3: Validate biological mechanism of detected variants.
    
    Returns: Dictionary with mechanism summary per gene
    """
    mechanism = {}
    
    for gene_name, gene_data in detected_genes.items():
        effects = []
        
        # Get gene info from database
        db_info = PIGMENTATION_GENES_DB.get(gene_name)
        
        if gene_data.high_impact > 0:
            # Analyze HIGH impact variants
            high_variants = [v for v in variants if v.gene == gene_name and v.impact == "HIGH"]
            consequences = set(v.consequence for v in high_variants)
            
            if "stop_gained" in consequences:
                effects.append(f"{sum(1 for v in high_variants if v.consequence == 'stop_gained')} stop_gained (truncated protein)")
            if "splice_donor_variant" in consequences or "splice_acceptor_variant" in consequences:
                splice_count = sum(1 for v in high_variants if "splice" in v.consequence)
                effects.append(f"{splice_count} splice variants (abnormal splicing)")
            if "frameshift_variant" in consequences:
                effects.append("frameshift (altered reading frame)")
        
        if gene_data.moderate_impact > 0:
            # Analyze MODERATE impact variants
            mod_variants = [v for v in variants if v.gene == gene_name and v.impact == "MODERATE"]
            missense_count = sum(1 for v in mod_variants if "missense" in v.consequence)
            if missense_count > 0:
                effects.append(f"{missense_count} missense variants (amino acid changes)")
        
        # Add functional context
        if db_info:
            function_summary = f"Function: {db_info.function[:80]}..."
        else:
            function_summary = "Unknown gene function"
        
        if effects:
            mechanism[gene_name] = f"{'; '.join(effects)}. {function_summary}"
        else:
            # Only LOW/MODIFIER variants
            mechanism[gene_name] = f"Regulatory/intronic variants only. {function_summary}"
    
    return mechanism


# ============================================================================
# GNOMAD API (OPTIONAL)
# ============================================================================

def fetch_gnomad_frequencies(rsids: List[str], timeout: int = 30) -> Dict[str, Dict]:
    """
    Fetch population allele frequencies from gnomAD API.
    
    Note: gnomAD API can be slow/rate-limited. This is optional.
    """
    if not HAS_REQUESTS:
        print("[WARN] requests module not available, skipping gnomAD")
        return {}
    
    frequencies = {}
    base_url = "https://gnomad.broadinstitute.org/api"
    
    # gnomAD GraphQL query
    query_template = """
    query getVariant($rsid: String!) {
        variant(rsid: $rsid, dataset: gnomad_r3) {
            rsid
            genome {
                populations {
                    id
                    ac
                    an
                    af
                }
            }
        }
    }
    """
    
    for rsid in rsids[:10]:  # Limit to 10 to avoid rate limiting
        try:
            response = requests.post(
                base_url,
                json={"query": query_template, "variables": {"rsid": rsid}},
                timeout=timeout,
                headers={"Content-Type": "application/json"}
            )
            
            if response.status_code == 200:
                data = response.json()
                variant_data = data.get("data", {}).get("variant")
                if variant_data and variant_data.get("genome"):
                    pops = variant_data["genome"].get("populations", [])
                    pop_freqs = {}
                    for pop in pops:
                        pop_id = pop.get("id", "")
                        af = pop.get("af")
                        if af is not None and pop_id in ["afr", "eur", "eas", "amr", "sas"]:
                            pop_freqs[pop_id.upper()] = af
                    if pop_freqs:
                        frequencies[rsid] = pop_freqs
            
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            print(f"[WARN] gnomAD fetch failed for {rsid}: {e}")
    
    return frequencies


# ============================================================================
# REPORT GENERATION
# ============================================================================

def calculate_validation_score(result: ValidationResult) -> Tuple[str, float]:
    """
    Calculate overall validation score and status.
    
    Scoring:
    - Gene match: 40% weight
    - Known rsIDs found: 20% weight  
    - Crawford overlap: 20% weight
    - Mechanism consistency: 20% weight
    
    Returns: (status, score 0-100)
    """
    score = 0.0
    
    # Gene match (40%)
    gene_score = result.gene_match_percentage * 0.4
    score += gene_score
    
    # Known rsIDs (20%) - any known rsID gives full points
    rsid_score = 20.0 if result.known_pigmentation_rsids_found else 0.0
    score += rsid_score
    
    # Crawford overlap (20%)
    crawford_score = result.crawford_overlap_percentage * 0.2
    score += crawford_score
    
    # Mechanism (20%) - based on functional variants found
    if result.total_high_impact > 0 and result.total_moderate_impact > 0:
        mechanism_score = 20.0
    elif result.total_high_impact > 0 or result.total_moderate_impact > 0:
        mechanism_score = 15.0
    elif result.total_missense > 0:
        mechanism_score = 10.0
    else:
        mechanism_score = 5.0
    score += mechanism_score
    
    # Determine status
    if score >= 80:
        status = "STRONGLY VALIDATED (+++)"
    elif score >= 60:
        status = "VALIDATED (++)"
    elif score >= 40:
        status = "PARTIALLY VALIDATED (+)"
    else:
        status = "NOT VALIDATED (-)"
    
    return status, score


def generate_report(result: ValidationResult, output_path: Path, detected_genes_data: Dict[str, DetectedGene]) -> None:
    """Generate the validation report in Markdown format."""
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Pigmentation Hypothesis Validation Report\n\n")
        f.write(f"> Generated: {timestamp}\n\n")
        
        # Executive Summary
        f.write("## Executive Summary\n\n")
        f.write("This report validates whether the DeepLIFT-identified genomic regions\n")
        f.write("correspond to known pigmentation genes, confirming the hypothesis that\n")
        f.write("the pipeline can identify genetic determinants of skin pigmentation.\n\n")
        f.write(f"### Validation Status: {result.validation_status}\n\n")
        
        # Validation Scores
        f.write("## Validation Scores\n\n")
        f.write("| Metric | Value |\n")
        f.write("|--------|-------|\n")
        f.write(f"| Total genes detected | {len(result.detected_genes)} |\n")
        f.write(f"| Primary pigmentation genes | {result.primary_genes_count} ({result.primary_genes_count/len(result.detected_genes)*100:.1f}%) |\n")
        all_pigment = result.primary_genes_count + result.regulatory_genes_count
        f.write(f"| All pigmentation genes (primary + regulatory) | {all_pigment} ({all_pigment/len(result.detected_genes)*100:.1f}%) |\n")
        extended = all_pigment + result.uv_response_genes_count
        f.write(f"| Extended (+ UV response) | {extended} ({extended/len(result.detected_genes)*100:.1f}%) |\n")
        f.write(f"| Unknown/unclassified | {len(result.unmatched_genes)} ({len(result.unmatched_genes)/len(result.detected_genes)*100:.1f}%) |\n\n")
        
        # Gene Classifications Table
        f.write("## Gene Classifications\n\n")
        f.write("| Gene | Category | OMIM | Function | HIGH Impact | MODERATE Impact |\n")
        f.write("|------|----------|------|----------|-------------|------------------|\n")
        
        for gene_name in result.detected_genes:
            gene_data = detected_genes_data.get(gene_name)
            db_info = PIGMENTATION_GENES_DB.get(gene_name)
            
            if db_info:
                category = db_info.category.replace("_", " ").title()
                omim = db_info.omim
                func = db_info.function[:50] + "..." if len(db_info.function) > 50 else db_info.function
            else:
                category = "Unknown"
                omim = "N/A"
                func = "Unknown function"
            
            high = gene_data.high_impact if gene_data else 0
            moderate = gene_data.moderate_impact if gene_data else 0
            
            f.write(f"| {gene_name} | {category} | {omim} | {func} | {high} | {moderate} |\n")
        
        f.write("\n")
        
        # Detailed Gene Analysis
        f.write("## Detailed Gene Analysis\n\n")
        
        for gene_name in result.detected_genes:
            gene_data = detected_genes_data.get(gene_name)
            db_info = PIGMENTATION_GENES_DB.get(gene_name)
            
            f.write(f"### {gene_name}\n\n")
            
            if db_info:
                f.write(f"**Category**: {db_info.category.replace('_', ' ').title()} Pigmentation Gene\n\n")
                f.write(f"**Function**: {db_info.function}\n\n")
                f.write(f"**OMIM**: [{db_info.omim}](https://omim.org/entry/{db_info.omim})\n\n")
                f.write(f"**GWAS Studies**: {', '.join(db_info.gwas_references)}\n\n")
                f.write(f"**Population Effect**: {db_info.population_effect}\n\n")
            else:
                f.write("**Category**: Unknown (not in pigmentation gene database)\n\n")
            
            if gene_data:
                f.write("**Variant Statistics**:\n")
                f.write(f"- Unique variants: {gene_data.unique_variants}\n")
                f.write(f"- HIGH impact: {gene_data.high_impact}\n")
                f.write(f"- MODERATE impact: {gene_data.moderate_impact}\n")
                f.write(f"- LOW impact: {gene_data.low_impact}\n")
                f.write(f"- Top consequences: {gene_data.top_consequences}\n\n")
        
        # Comparison with Published GWAS
        f.write("## Comparison with Published GWAS\n\n")
        f.write("### Crawford et al., 2017 (Science)\n\n")
        f.write("This landmark study identified key pigmentation genes in African populations:\n\n")
        f.write("| Crawford et al. Gene | Detected by Pipeline |\n")
        f.write("|---------------------|---------------------|\n")
        detected_set = set(result.detected_genes)
        for gene in sorted(CRAWFORD_2017_GENES):
            status = "Yes" if gene in detected_set else "No"
            f.write(f"| {gene} | {status} |\n")
        f.write(f"\n**Overlap with Crawford et al.**: {result.crawford_overlap_count}/{len(CRAWFORD_2017_GENES)} ({result.crawford_overlap_percentage:.1f}%)\n\n")
        
        # Known rsIDs found
        if result.known_pigmentation_rsids_found:
            f.write("## Known Pigmentation rsIDs Detected\n\n")
            f.write("The following well-characterized pigmentation SNPs were found in the data:\n\n")
            f.write("| rsID | Gene | Effect | Reference |\n")
            f.write("|------|------|--------|----------|\n")
            for rsid in result.known_pigmentation_rsids_found:
                info = KNOWN_PIGMENTATION_RSIDS.get(rsid, {})
                gene = info.get("gene", "Unknown")
                effect = info.get("effect", "Unknown")[:50] + "..."
                refs = ", ".join(info.get("references", ["N/A"]))
                f.write(f"| {rsid} | {gene} | {effect} | {refs} |\n")
            f.write("\n")
        
        # Mechanism Summary
        if result.mechanism_summary:
            f.write("## Biological Mechanism Analysis\n\n")
            for gene, mechanism in result.mechanism_summary.items():
                f.write(f"**{gene}**: {mechanism}\n\n")
        
        # Population Frequencies (if available)
        if result.population_frequencies:
            f.write("## Population Allele Frequencies (gnomAD)\n\n")
            f.write("| rsID | AFR | EUR | EAS | AMR | SAS |\n")
            f.write("|------|-----|-----|-----|-----|-----|\n")
            for rsid, freqs in result.population_frequencies.items():
                afr = freqs.get("AFR", "N/A")
                eur = freqs.get("EUR", "N/A")
                eas = freqs.get("EAS", "N/A")
                amr = freqs.get("AMR", "N/A")
                sas = freqs.get("SAS", "N/A")
                
                afr_str = f"{afr:.4f}" if isinstance(afr, float) else afr
                eur_str = f"{eur:.4f}" if isinstance(eur, float) else eur
                eas_str = f"{eas:.4f}" if isinstance(eas, float) else eas
                amr_str = f"{amr:.4f}" if isinstance(amr, float) else amr
                sas_str = f"{sas:.4f}" if isinstance(sas, float) else sas
                
                f.write(f"| {rsid} | {afr_str} | {eur_str} | {eas_str} | {amr_str} | {sas_str} |\n")
            f.write("\n")
        
        # Conclusion
        f.write("## Conclusion\n\n")
        if "STRONGLY" in result.validation_status:
            f.write("The hypothesis is **strongly validated**. The DeepLIFT pipeline successfully\n")
            f.write("identified genes that are well-established in the pigmentation literature.\n")
            f.write("This provides strong evidence that the same approach can be applied to\n")
            f.write("identify longevity-associated genes when appropriate training data is available.\n")
        elif "VALIDATED" in result.validation_status and "NOT" not in result.validation_status:
            f.write("The hypothesis is **validated**. The pipeline identified relevant pigmentation genes,\n")
            f.write("supporting its potential application for other phenotype studies.\n")
        elif "PARTIALLY" in result.validation_status:
            f.write("The hypothesis is **partially validated**. Some pigmentation genes were identified,\n")
            f.write("but further optimization may improve detection.\n")
        else:
            f.write("The hypothesis **was not validated**. The detected genes do not correspond\n")
            f.write("to known pigmentation genes. Review pipeline parameters and input data.\n")
        
        f.write("\n---\n\n")
        f.write("*Report generated by validate_pigmentation_hypothesis.py*")
    
    print(f"[SAVED] Validation report: {output_path}")


# ============================================================================
# MAIN
# ============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate pigmentation hypothesis using DeepLIFT+VEP results"
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing annotate_deeplift_windows.py output files"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Output report file (default: <input_dir>/pigmentation_validation.md)"
    )
    parser.add_argument(
        "--gnomad",
        action="store_true",
        help="Fetch population frequencies from gnomAD API (optional, slow)"
    )
    
    args = parser.parse_args()
    
    input_dir = args.input_dir
    if not input_dir.is_dir():
        print(f"[ERROR] Input directory not found: {input_dir}", file=sys.stderr)
        return 1
    
    # Define file paths
    summary_by_gene_path = input_dir / "summary_by_gene.tsv"
    vep_variants_path = input_dir / "vep_annotated_variants.tsv"
    
    # Check required files
    if not summary_by_gene_path.exists():
        print(f"[ERROR] Required file not found: {summary_by_gene_path}", file=sys.stderr)
        return 1
    
    print(f"[INFO] Reading data from: {input_dir}")
    
    # Parse input files
    detected_genes = parse_summary_by_gene(summary_by_gene_path)
    print(f"[INFO] Detected genes: {', '.join(detected_genes.keys())}")
    
    variants = []
    if vep_variants_path.exists():
        variants = parse_vep_annotated_variants(vep_variants_path)
        print(f"[INFO] Loaded {len(variants)} VEP-annotated variants")
    else:
        print(f"[WARN] VEP variants file not found: {vep_variants_path}")
    
    # Level 1: Gene validation
    print("\n=== Level 1: Gene Validation ===")
    matched, unmatched, primary, regulatory, uv, crawford_count, crawford_pct = validate_genes(detected_genes)
    print(f"  Matched pigmentation genes: {matched}")
    print(f"  Unmatched genes: {unmatched}")
    print(f"  Crawford et al. 2017 overlap: {crawford_count}/{len(CRAWFORD_2017_GENES)} ({crawford_pct:.1f}%)")
    
    gene_match_pct = len(matched) / len(detected_genes) * 100 if detected_genes else 0
    
    # Level 2: Variant validation
    print("\n=== Level 2: Variant Validation ===")
    all_rsids, known_rsids, high_count, mod_count, missense_count = validate_variants(detected_genes, variants)
    print(f"  Total rsIDs found: {len(all_rsids)}")
    print(f"  Known pigmentation rsIDs: {known_rsids if known_rsids else 'None'}")
    print(f"  HIGH impact variants: {high_count}")
    print(f"  MODERATE impact variants: {mod_count}")
    print(f"  Missense variants: {missense_count}")
    
    # Level 3: Mechanism validation
    print("\n=== Level 3: Mechanism Validation ===")
    mechanism = validate_mechanism(detected_genes, variants)
    for gene, summary in mechanism.items():
        print(f"  {gene}: {summary[:80]}...")
    
    # Optional: gnomAD frequencies
    pop_frequencies = None
    if args.gnomad and all_rsids:
        print("\n=== Fetching gnomAD Frequencies ===")
        pop_frequencies = fetch_gnomad_frequencies(all_rsids)
        if pop_frequencies:
            print(f"  Fetched frequencies for {len(pop_frequencies)} variants")
        else:
            print("  No frequencies retrieved")
    
    # Create result object
    result = ValidationResult(
        detected_genes=list(detected_genes.keys()),
        matched_pigmentation_genes=matched,
        unmatched_genes=unmatched,
        gene_match_percentage=gene_match_pct,
        primary_genes_count=primary,
        regulatory_genes_count=regulatory,
        uv_response_genes_count=uv,
        crawford_overlap_count=crawford_count,
        crawford_overlap_percentage=crawford_pct,
        all_rsids_found=all_rsids,
        known_pigmentation_rsids_found=known_rsids,
        total_high_impact=high_count,
        total_moderate_impact=mod_count,
        total_missense=missense_count,
        mechanism_summary=mechanism,
        population_frequencies=pop_frequencies
    )
    
    # Calculate validation score
    status, score = calculate_validation_score(result)
    result.validation_status = status
    result.validation_score = score
    
    print(f"\n=== Validation Result ===")
    print(f"  Status: {status}")
    print(f"  Score: {score:.1f}/100")
    
    # Generate report
    output_path = args.output or (input_dir / "pigmentation_validation.md")
    generate_report(result, output_path, detected_genes)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

