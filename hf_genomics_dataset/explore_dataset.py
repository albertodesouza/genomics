#!/usr/bin/env python3
"""
explore_dataset.py

Script rápido para explorar e filtrar o dataset HuggingFace Arrow
gerado pelo build_dataset.py, com suporte a views JSON.

Uso sem view (exibe tudo):
    python3 hf_genomics_dataset/explore_dataset.py

Uso com view:
    python3 hf_genomics_dataset/explore_dataset.py --view hf_genomics_dataset/views/mc1r_afr.view.json
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from datasets import load_from_disk


# ── Argumentos ───────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Explorar dataset HF com view opcional.")
parser.add_argument("--dataset", default="/home/felipe/final_test", help="Diretório do dataset Arrow.")
parser.add_argument("--view", default=None, help="Caminho para um arquivo .view.json.")
args = parser.parse_args()

OUTPUT_DIR = Path(args.dataset)

# ── Carrega view (se fornecida) ──────────────────────────────────────────────
view = {}
if args.view:
    view_path = Path(args.view)
    if not view_path.exists():
        raise FileNotFoundError(f"View não encontrada: {view_path}")
    with open(view_path) as f:
        view = json.load(f)
    print(f"=== View carregada: '{view.get('name', view_path.stem)}' ===")
    print(f"    {view.get('description', '')}")
    # Se a view definir dataset_dir, usa ele
    if "dataset_dir" in view:
        OUTPUT_DIR = Path(view["dataset_dir"])

genes_to_use        = view.get("genes_to_use")          # ex: ["MC1R", "TYR"]
superpops_to_use    = set(view.get("superpopulations_to_use") or [])
pops_to_use         = set(view.get("populations_to_use") or [])
sample_ids_to_use   = set(view.get("sample_ids") or [])

# ── Carrega as três tabelas ──────────────────────────────────────────────────
print(f"\n=== Carregando dataset de {OUTPUT_DIR} ===")
samples        = load_from_disk(str(OUTPUT_DIR / "samples"))
gene_windows   = load_from_disk(str(OUTPUT_DIR / "gene_windows"))
gene_reference = load_from_disk(str(OUTPUT_DIR / "gene_reference"))

print(f"Samples (total):       {len(samples)} linhas")
print(f"Gene windows (total):  {len(gene_windows)} linhas")
print(f"Gene references (total): {len(gene_reference)} linhas")

# ── Aplica filtros da view ───────────────────────────────────────────────────
if sample_ids_to_use:
    samples      = samples.filter(lambda x: x["sample_id"] in sample_ids_to_use)
    gene_windows = gene_windows.filter(lambda x: x["sample_id"] in sample_ids_to_use)

if superpops_to_use:
    samples      = samples.filter(lambda x: x["superpopulation"] in superpops_to_use)
    gene_windows = gene_windows.filter(lambda x: x["sample_id"] in set(samples["sample_id"]))

if pops_to_use:
    samples      = samples.filter(lambda x: x["population"] in pops_to_use)
    gene_windows = gene_windows.filter(lambda x: x["sample_id"] in set(samples["sample_id"]))

if genes_to_use:
    gene_windows = gene_windows.filter(lambda x: x["gene"] in genes_to_use)
    gene_reference = gene_reference.filter(lambda x: x["gene_symbol"] in genes_to_use)

print(f"\nSamples (após filtro):       {len(samples)} linhas")
print(f"Gene windows (após filtro):  {len(gene_windows)} linhas")
print(f"Gene reference (após filtro): {len(gene_reference)} linhas")

# ── Inspeciona uma amostra ───────────────────────────────────────────────────
if len(samples) > 0:
    print("\n=== Primeira linha de samples ===")
    for k, v in samples[0].items():
        print(f"  {k}: {str(v)[:120]}")

# Build a quick lookup for references
ref_lookup = {row["gene_symbol"]: row for row in gene_reference}

import io
import numpy as np

def format_value(v):
    if isinstance(v, bytes):
        size_kb = len(v) / 1024
        size_str = f"{size_kb/1024:.2f} MB" if size_kb > 1024 else f"{size_kb:.2f} KB"
        
        # Tenta decodificar se for maior que 0
        if len(v) > 0:
            try:
                data = np.load(io.BytesIO(v))
                if "values" in data:
                    arr = data["values"]
                    # Pega um pequeno slice central ou inicial para mostrar
                    preview = arr.flatten()[:5]
                    return f"<binary: {size_str}, shape: {arr.shape}, preview: {preview}...>"
            except Exception:
                pass
        return f"<binary: {size_str}>"
    if isinstance(v, dict):
        return {k: format_value(v2) for k, v2 in v.items()}
    return str(v)[:120]

if len(gene_windows) > 0:
    for i in range(0, min(11, len(gene_windows))):
        print(f"\n=== Linha {i} de gene_windows ===")
        window = gene_windows[i]
        for k, v in window.items():
            print(f"  {k}: {format_value(v)}")
        
        gene_symbol = window["gene"]
        if gene_symbol in ref_lookup:
            ref = ref_lookup[gene_symbol]
            print(f"  [REFERENCE TABLE JOIN]")
            print(f"    ref_sequence preview: {ref['ref_sequence'][:80]}...")

# ── Estatísticas ─────────────────────────────────────────────────────────────
print("\n=== Genes disponíveis (após filtro) ===")
print(f"  {sorted(set(gene_windows['gene']))}")

print("\n=== Distribuição por superpopulação ===")
for pop, count in sorted(Counter(samples["superpopulation"]).items()):
    print(f"  {pop}: {count}")

# ── Prévia como DataFrame ────────────────────────────────────────────────────
print("\n=== Prévia como DataFrame ===")
df = samples.to_pandas()[["sample_id", "population", "superpopulation", "longevity", "available_genes"]]
print(df.to_string(index=False))

