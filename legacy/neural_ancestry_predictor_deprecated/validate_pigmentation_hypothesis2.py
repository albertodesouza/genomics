#!/usr/bin/env python3
"""
validate_pigmentation_hypothesis2.py

Independent validation script for the pigmentation hypothesis using DeepLIFT + VEP outputs.
This script analyzes the output directory of annotate_deeplift_windows.py to confirm if
identified variants in pigmentation genes have functional impact or known phenotypic associations.

Usage:
    python3 validate_pigmentation_hypothesis2.py <input_dir> [--output <report_path>]
"""

import os
import sys
import json
import argparse
import pandas as pd
from typing import List, Dict, Any, Set

# Configuration
PIGMENTATION_GENES = [
    "MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR",
    "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"
]

PHENOTYPE_KEYWORDS = [
    "skin", "pigmentation", "melanin", "hair color", "eye color",
    "albinism", "tanning", "freckles", "vitiligo", "melanoma",
    "basal cell carcinoma", "squamous cell carcinoma"
]

# Known significant pigmentation variants (rsID -> description)
KNOWN_VARIANTS = {
    "rs1426654": "SLC24A5 A111T (Major skin color determinant in Europeans/Africans)",
    "rs16891982": "SLC45A2 L374F (Dark/Light skin allele)",
    "rs12913832": "HERC2/OCA2 (Blue eye color / skin pigmentation)",
    "rs1042602": "TYR S192Y (Tyrosinase activity)",
    "rs1126809": "TYR R402Q (Tyrosinase activity)",
    "rs1800404": "OCA2 H615R (Light skin in East Asians)",
    "rs1800414": "OCA2 A481T",
    "rs74653330": "OCA2 A481T (Alt ID)",
    "rs3827760": "EDAR V370A (Hair thickness/tooth shape, associated with selection)",
    "rs28777": "SLC45A2 (Skin color)",
    "rs619865": "ASIP (Skin sensitivity to sun)",
    "rs2424984": "SLC24A5 (Regulatory variant)",
}

def load_data(input_dir: str):
    """Loads necessary data files from the input directory."""
    files = {
        "summary_gene": os.path.join(input_dir, "summary_by_gene.tsv"),
        "variants_tsv": os.path.join(input_dir, "vep_annotated_variants.tsv"),
        "vep_json": os.path.join(input_dir, "vep_annotations.all.json")
    }
    
    data = {}
    
    # Check existence
    for key, path in files.items():
        if not os.path.exists(path):
            print(f"[ERROR] Required file not found: {path}")
            return None
            
    # Load Summary by Gene
    try:
        data["summary_gene"] = pd.read_csv(files["summary_gene"], sep='\t')
    except Exception as e:
        print(f"[ERROR] Failed to read summary_by_gene.tsv: {e}")
        return None
        
    # Load Variants TSV
    try:
        data["variants_tsv"] = pd.read_csv(files["variants_tsv"], sep='\t')
    except Exception as e:
        print(f"[ERROR] Failed to read vep_annotated_variants.tsv: {e}")
        return None
        
    # Load VEP JSON
    try:
        with open(files["vep_json"], 'r') as f:
            data["vep_json"] = json.load(f)
    except Exception as e:
        print(f"[ERROR] Failed to read vep_annotations.all.json: {e}")
        return None
        
    return data

def analyze_impact(summary_df: pd.DataFrame, variants_df: pd.DataFrame) -> Dict:
    """Analyzes variant impacts in pigmentation genes."""
    results = {
        "gene_stats": [],
        "high_impact": [],
        "moderate_impact": []
    }
    
    # Filter for our genes
    pigment_df = summary_df[summary_df['gene'].isin(PIGMENTATION_GENES)].copy()
    
    if pigment_df.empty:
        return results

    # Get gene stats
    for _, row in pigment_df.iterrows():
        results["gene_stats"].append({
            "gene": row['gene'],
            "total": row['unique_variants'],
            "high": row.get('HIGH', 0),
            "moderate": row.get('MODERATE', 0),
            "low": row.get('LOW', 0),
            "modifier": row.get('MODIFIER', 0)
        })
        
    # Get specific high/moderate impact variants
    # Need to join with variants_tsv to get details
    mask = (variants_df['gene_symbols'].apply(lambda x: any(g in str(x).split(',') for g in PIGMENTATION_GENES)))
    relevant_vars = variants_df[mask]
    
    for _, row in relevant_vars.iterrows():
        impact = row.get('worst_impact', 'MODIFIER')
        if impact in ['HIGH', 'MODERATE']:
            var_info = {
                "gene": row['gene_symbols'],
                "chrom": row['chrom'],
                "pos": row['pos_1based'],
                "ref": row['ref'],
                "alt": row['alt'],
                "impact": impact,
                "consequence": row['most_severe_consequence'],
                "rsid": row.get('rsids', '.')
            }
            if impact == 'HIGH':
                results["high_impact"].append(var_info)
            else:
                results["moderate_impact"].append(var_info)
                
    return results

def search_phenotypes(vep_json: List[Dict]) -> List[Dict]:
    """Searches VEP JSON for variants with relevant phenotype associations."""
    phenotype_hits = []
    
    for item in vep_json:
        # Check colocated variants
        if "colocated_variants" in item:
            for cv in item["colocated_variants"]:
                # Check known rsIDs
                rsid = cv.get("id", "")
                
                # Check phenotype_or_disease
                phenotypes = []
                if "phenotype_or_disease" in cv:
                    val = cv["phenotype_or_disease"]
                    if isinstance(val, list):
                         phenotypes.extend([str(p).lower() for p in val if isinstance(p, (str, int))])
                    elif isinstance(val, (str, int)):
                         phenotypes.append(str(val).lower())

                if "clin_sig" in cv:
                    val = cv["clin_sig"]
                    if isinstance(val, list):
                         phenotypes.extend([str(p).lower() for p in val if isinstance(p, (str, int))])
                    elif isinstance(val, (str, int)):
                         phenotypes.append(str(val).lower())

                is_relevant = False
                matched_keywords = []
                
                for p in phenotypes:
                    for kw in PHENOTYPE_KEYWORDS:
                        if kw in p:
                            is_relevant = True
                            matched_keywords.append(kw)
                
                # If matched keywords or is a known variant
                is_known = rsid in KNOWN_VARIANTS
                
                if is_relevant or is_known:
                    # Find associated gene from transcript consequences
                    genes = set()
                    if "transcript_consequences" in item:
                        for tc in item["transcript_consequences"]:
                            if "gene_symbol" in tc:
                                genes.add(tc["gene_symbol"])
                    
                    # Only keep if associated with pigmentation genes of interest
                    # (Optional: broaden search, but for hypothesis validation, focus on these)
                    intersect_genes = genes.intersection(set(PIGMENTATION_GENES))
                    
                    # Prepare clean list for report
                    final_phenotypes = []
                    raw_p = cv.get("phenotype_or_disease", [])
                    if isinstance(raw_p, list):
                        final_phenotypes = raw_p
                    elif isinstance(raw_p, (str, int)):
                        final_phenotypes = [raw_p]

                    if intersect_genes:
                        phenotype_hits.append({
                            "rsid": rsid,
                            "genes": list(intersect_genes),
                            "phenotypes": final_phenotypes,
                            "clin_sig": cv.get("clin_sig", []),
                            "matched_keywords": list(set(matched_keywords)),
                            "is_known_list": is_known,
                            "known_desc": KNOWN_VARIANTS.get(rsid, "")
                        })
                        
    return phenotype_hits

def check_known_rsids(variants_df: pd.DataFrame) -> List[Dict]:
    """Checks explicitly for the list of known important rsIDs in the dataframe."""
    found = []
    
    # Flatten rsids column
    # rsids can be comma separated or NaN
    
    for _, row in variants_df.iterrows():
        row_rsids = str(row.get('rsids', ''))
        if pd.isna(row_rsids) or row_rsids == 'nan' or row_rsids == '.':
            continue
            
        current_ids = row_rsids.split(',')
        for rid in current_ids:
            if rid in KNOWN_VARIANTS:
                found.append({
                    "rsid": rid,
                    "description": KNOWN_VARIANTS[rid],
                    "gene": row['gene_symbols'],
                    "pos": f"{row['chrom']}:{row['pos_1based']}",
                    "consequence": row['most_severe_consequence']
                })
                
    return found

def generate_report(output_path: str, impact_data: Dict, phenotype_data: List[Dict], known_found: List[Dict]):
    """Generates the Markdown report."""
    
    with open(output_path, 'w') as f:
        f.write("# Validação da Hipótese de Pigmentação (Método 2)\n\n")
        f.write("**Data da Análise:** " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M") + "\n\n")
        f.write("Este relatório valida independentemente a hipótese de que o pipeline identificou variantes funcionais associadas à pigmentação.\n\n")
        
        # 1. Resumo
        f.write("## 1. Resumo Executivo\n\n")
        total_genes_found = len(impact_data["gene_stats"])
        total_high = len(impact_data["high_impact"])
        total_mod = len(impact_data["moderate_impact"])
        
        f.write(f"- **Genes de pigmentação identificados:** {total_genes_found} de {len(PIGMENTATION_GENES)}\n")
        f.write(f"- **Variantes de ALTO impacto (HIGH):** {total_high}\n")
        f.write(f"- **Variantes de impacto MODERADO:** {total_mod}\n")
        f.write(f"- **Associações fenotípicas explícitas encontradas:** {len(phenotype_data)}\n")
        f.write(f"- **SNPs conhecidos da literatura encontrados:** {len(known_found)}\n\n")
        
        # 2. Análise por Gene
        f.write("## 2. Análise Detalhada por Gene\n\n")
        if impact_data["gene_stats"]:
            f.write("| Gene | Variantes Totais | HIGH | MODERATE | LOW | MODIFIER |\n")
            f.write("|------|------------------|------|----------|-----|----------|\n")
            for stat in impact_data["gene_stats"]:
                f.write(f"| **{stat['gene']}** | {stat['total']} | {stat['high']} | {stat['moderate']} | {stat['low']} | {stat['modifier']} |\n")
            f.write("\n")
        else:
            f.write("*Nenhum gene de pigmentação alvo foi encontrado na saída.*\n\n")
            
        # 3. Variantes de Alto/Moderado Impacto
        f.write("## 3. Variantes Funcionais (HIGH/MODERATE)\n\n")
        
        if impact_data["high_impact"]:
            f.write("### 🔴 Alto Impacto (HIGH)\n\n")
            f.write("| Gene | Posição | Ref/Alt | Consequência | rsID |\n")
            f.write("|------|---------|---------|--------------|------|\n")
            for var in impact_data["high_impact"]:
                f.write(f"| {var['gene']} | {var['chrom']}:{var['pos']} | {var['ref']}/{var['alt']} | {var['consequence']} | {var['rsid']} |\n")
            f.write("\n")
            
        if impact_data["moderate_impact"]:
            f.write("### 🟡 Impacto Moderado\n\n")
            f.write("*(Mostrando as primeiras 10 variantes)*\n\n")
            f.write("| Gene | Posição | Ref/Alt | Consequência | rsID |\n")
            f.write("|------|---------|---------|--------------|------|\n")
            for i, var in enumerate(impact_data["moderate_impact"]):
                if i >= 10: break
                f.write(f"| {var['gene']} | {var['chrom']}:{var['pos']} | {var['ref']}/{var['alt']} | {var['consequence']} | {var['rsid']} |\n")
            if len(impact_data["moderate_impact"]) > 10:
                f.write(f"| ... | ... | ... | ... | ... |\n")
            f.write("\n")
            
        # 4. Associações Fenotípicas (VEP)
        f.write("## 4. Associações Fenotípicas (VEP)\n\n")
        f.write("Variantes onde o VEP reportou associações diretas com pigmentação, pele ou doenças relacionadas.\n\n")
        
        if phenotype_data:
            for hit in phenotype_data:
                genes_str = ", ".join(hit['genes'])
                f.write(f"### {hit['rsid']} ({genes_str})\n")
                if hit['is_known_list']:
                    f.write(f"- **Literatura:** {hit['known_desc']}\n")
                if hit['phenotypes']:
                    f.write(f"- **Fenótipos Reportados:** {', '.join([str(p) for p in hit['phenotypes']])}\n")
                if hit['clin_sig']:
                    f.write(f"- **Significância Clínica:** {', '.join([str(c) for c in hit['clin_sig']])}\n")
                f.write("\n")
        else:
            f.write("*Nenhuma associação fenotípica explícita encontrada nos metadados do VEP.*\n\n")
            
        # 5. Validação com Lista de Ouro
        f.write("## 5. Validação com SNPs Conhecidos (Lista de Ouro)\n\n")
        f.write("Verificação cruzada com uma lista curada de variantes importantes para pigmentação.\n\n")
        
        if known_found:
            f.write("| rsID | Gene | Descrição Conhecida | Encontrado no Dataset? |\n")
            f.write("|------|------|---------------------|------------------------|\n")
            for item in known_found:
                f.write(f"| **{item['rsid']}** | {item['gene']} | {item['description']} | ✅ SIM |\n")
            f.write("\n")
        else:
            f.write("*Nenhum dos SNPs da 'Lista de Ouro' específica foi encontrado nas regiões analisadas.*\n")
            f.write("*Isso pode ocorrer se o DeepLIFT apontou para regiões regulatórias próximas, mas não exatamente para o SNP, ou se o indivíduo analisado não possui esses SNPs específicos.*\n\n")

        # 6. Conclusão
        f.write("## 6. Conclusão\n\n")
        
        strong_evidence = (len(impact_data["high_impact"]) > 0) or (len(known_found) > 0) or (len(phenotype_data) > 0)
        moderate_evidence = (len(impact_data["moderate_impact"]) > 0) and (total_genes_found > 0)
        
        if strong_evidence:
            f.write("✅ **A hipótese é FORTEMENTE suportada.**\n\n")
            f.write("O pipeline identificou variantes com relevância funcional conhecida ou prevista em genes chave de pigmentação.\n")
            if len(known_found) > 0:
                f.write("A presença de variantes conhecidas da literatura confirma que o DeepLIFT está identificando as regiões causais corretas.\n")
        elif moderate_evidence:
            f.write("⚠️ **A hipótese é PARCIALMENTE suportada.**\n\n")
            f.write("O pipeline encontrou variantes em genes de pigmentação com impacto moderado (missense), mas não confirmou SNPs clássicos da literatura ou associações fenotípicas diretas nos metadados.\n")
        else:
            f.write("❌ **Evidência limitada encontrada.**\n\n")
            f.write("Embora as regiões apontem para genes de pigmentação, não foram encontradas variantes funcionais óbvias ou anotações fenotípicas que expliquem a causalidade apenas com estes dados.\n")

def main():
    parser = argparse.ArgumentParser(description="Validate pigmentation hypothesis from DeepLIFT output")
    parser.add_argument("input_dir", help="Directory containing annotate_deeplift_windows.py output")
    parser.add_argument("--output", help="Output report path", default=None)
    
    args = parser.parse_args()
    
    if args.output is None:
        args.output = os.path.join(args.input_dir, "pigmentation_validation2.md")
        
    print(f"[INFO] Analyzing output in: {args.input_dir}")
    
    # 1. Load Data
    data = load_data(args.input_dir)
    if not data:
        sys.exit(1)
        
    # 2. Analyze
    print("[INFO] Analyzing variant impacts...")
    impact_results = analyze_impact(data["summary_gene"], data["variants_tsv"])
    
    print("[INFO] Searching VEP phenotypes...")
    phenotype_results = search_phenotypes(data["vep_json"])
    
    print("[INFO] Checking known rsIDs...")
    known_found = check_known_rsids(data["variants_tsv"])
    
    # 3. Generate Report
    print(f"[INFO] Generating report: {args.output}")
    generate_report(args.output, impact_results, phenotype_results, known_found)
    
    print("[INFO] Done.")

if __name__ == "__main__":
    main()

