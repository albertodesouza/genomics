"""Well-characterised non-coding regulatory elements, offered as window presets next to genes.

A gene window is centred on the gene body, so a distal enhancer can sit near the edge of (or
outside) a 16-131 kb window. These presets centre the window on a causal or lead variant in the
element instead. Positions are GRCh38, 1-based, checked against Ensembl (release 113+).

Each entry becomes an import/prediction region ``{name, chrom, start, end}`` with ``start == end``
at the variant.
"""
from __future__ import annotations

from typing import Any, Dict, List

REGULATORY_REGIONS: List[Dict[str, Any]] = [
    # -------------------------------------------------------------- pigmentation (melanocytes)
    {
        "name": "OCA2_enh_rs12913832", "chrom": "chr15", "position": 28120472, "rsid": "rs12913832",
        "category": "Pigmentation", "targets": ["OCA2"],
        "description": "HERC2 intron 86 enhancer that loops to the OCA2 promoter; rs12913832 is the main determinant of blue vs brown eyes and modulates skin/hair pigmentation.",
        "reference": "Visser et al. 2012, Genome Res; Sturm et al. 2008, Am J Hum Genet",
    },
    {
        "name": "KITLG_enh_rs12821256", "chrom": "chr12", "position": 88934558, "rsid": "rs12821256",
        "category": "Pigmentation", "targets": ["KITLG"],
        "description": "Hair-follicle enhancer ~350 kb upstream of KITLG; rs12821256 weakens a LEF1 site and causes blond hair in Northern Europeans.",
        "reference": "Guenther et al. 2014, Nat Genet",
    },
    {
        "name": "IRF4_enh_rs12203592", "chrom": "chr6", "position": 396321, "rsid": "rs12203592",
        "category": "Pigmentation", "targets": ["IRF4", "TYR"],
        "description": "IRF4 intron 4 melanocyte enhancer; rs12203592 disrupts TFAP2A-dependent activation (freckling, hair colour, tanning).",
        "reference": "Praetorius et al. 2013, Cell",
    },
    {
        "name": "BNC2_enh_rs12350739", "chrom": "chr9", "position": 16885019, "rsid": "rs12350739",
        "category": "Pigmentation", "targets": ["BNC2"],
        "description": "Intergenic enhancer ~13 kb upstream of BNC2; rs12350739 changes its activity and BNC2 expression in melanocytes (skin colour saturation).",
        "reference": "Visser et al. 2014, Hum Mol Genet",
    },
    {
        "name": "TMEM138_DDB1_rs7948623", "chrom": "chr11", "position": 61369675, "rsid": "rs7948623",
        "category": "Pigmentation", "targets": ["TMEM138", "DDB1"],
        "description": "Regulatory element between TMEM138 and DDB1 associated with skin pigmentation in African populations; eQTL for both genes.",
        "reference": "Crawford et al. 2017, Science",
    },
    # ---------------------------------------------------------------- classic regulatory loci
    {
        "name": "LCT_enh_rs4988235", "chrom": "chr2", "position": 135851076, "rsid": "rs4988235",
        "category": "Classic enhancers", "targets": ["LCT"],
        "description": "MCM6 intron 13 enhancer of LCT; the -13910 C>T allele (rs4988235) confers lactase persistence in Europeans.",
        "reference": "Enattah et al. 2002, Nat Genet",
    },
    {
        "name": "FTO_IRX3_rs1421085", "chrom": "chr16", "position": 53767042, "rsid": "rs1421085",
        "category": "Classic enhancers", "targets": ["IRX3", "IRX5"],
        "description": "FTO intron 1 enhancer controlling IRX3/IRX5 in adipocyte progenitors; rs1421085 disrupts an ARID5B motif (obesity risk).",
        "reference": "Claussnitzer et al. 2015, N Engl J Med",
    },
    {
        "name": "MYC_8q24_rs6983267", "chrom": "chr8", "position": 127401060, "rsid": "rs6983267",
        "category": "Classic enhancers", "targets": ["MYC"],
        "description": "MYC-335 enhancer in the 8q24 gene desert, >300 kb from MYC; rs6983267 alters a TCF7L2 site (colorectal and prostate cancer risk).",
        "reference": "Pomerantz et al. 2009, Nat Genet; Tuupanen et al. 2009, Nat Genet",
    },
    {
        "name": "BCL11A_enh_rs1427407", "chrom": "chr2", "position": 60490908, "rsid": "rs1427407",
        "category": "Classic enhancers", "targets": ["BCL11A"],
        "description": "+62 site of the erythroid BCL11A intron 2 enhancer; rs1427407 sets fetal haemoglobin levels (the adjacent +58 site is the CRISPR therapy target).",
        "reference": "Bauer et al. 2013, Science",
    },
    {
        "name": "SORT1_enh_rs12740374", "chrom": "chr1", "position": 109274968, "rsid": "rs12740374",
        "category": "Classic enhancers", "targets": ["SORT1"],
        "description": "1p13 liver enhancer (CELSR2 3' end); rs12740374 creates a C/EBP site that raises SORT1 expression and lowers LDL cholesterol.",
        "reference": "Musunuru et al. 2010, Nature",
    },
    {
        "name": "CDKN2B_9p21_rs1333049", "chrom": "chr9", "position": 22125504, "rsid": "rs1333049",
        "category": "Classic enhancers", "targets": ["CDKN2A", "CDKN2B"],
        "description": "9p21 coronary artery disease interval (ANRIL/CDKN2B-AS1), dense in enhancers that interact with CDKN2A/B; rs1333049 is a lead risk SNP.",
        "reference": "WTCCC 2007, Nature; Harismendy et al. 2011, Nature",
    },
]


def preset_regions() -> List[Dict[str, Any]]:
    """Presets with the ``{name, chrom, start, end}`` keys the importer takes."""
    return [{**r, "start": r["position"], "end": r["position"]} for r in REGULATORY_REGIONS]
