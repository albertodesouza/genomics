#!/usr/bin/env python3
import subprocess as sp
from pathlib import Path

print("Baixando arquivo do 1000 Genomes para análise...")
url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"

temp_file = Path("temp_samples.tsv")
sp.run(['wget', '-q', '-O', str(temp_file), url], check=True)

print("\n📋 Primeiras 3 linhas do arquivo:")
print("=" * 80)

with open(temp_file, 'r') as f:
    for i, line in enumerate(f):
        if i >= 3:
            break
        print(f"Linha {i}:")
        fields = line.strip().split('\t')
        for j, field in enumerate(fields):
            if len(field) > 50:
                field = field[:47] + "..."
            print(f"  Col {j:2d}: {field}")
        print()

print("=" * 80)
print("\n🔍 Procurando por IDs de amostra (NA* ou HG*)...")

sample_ids = []
with open(temp_file, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:  # Pular primeira linha
            continue
        fields = line.strip().split('\t')
        for j, field in enumerate(fields):
            if field.startswith('NA') or field.startswith('HG'):
                sample_ids.append((j, field))
                if len(sample_ids) <= 5:
                    print(f"  Coluna {j}: {field}")
                if len(sample_ids) >= 10:
                    break
        if len(sample_ids) >= 10:
            break

if sample_ids:
    col_index = sample_ids[0][0]
    print(f"\n✓ IDs de amostra encontrados na coluna {col_index}")
else:
    print("\n✗ Nenhum ID encontrado")

temp_file.unlink()
print("\n✓ Arquivo temporário removido")
