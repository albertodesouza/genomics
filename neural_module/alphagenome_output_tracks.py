import os
import pandas as pd
from alphagenome.models import dna_client as alphagenome_client

# Configurar pandas para exibir todos os dados
pd.set_option('display.max_rows', None)       # Mostrar todas as linhas
pd.set_option('display.max_columns', None)    # Mostrar todas as colunas
pd.set_option('display.width', None)          # Sem limite de largura
pd.set_option('display.max_colwidth', None)   # Sem limite de largura de coluna

api_key = os.environ.get('ALPHAGENOME_API_KEY')
if not api_key:
    raise ValueError("ALPHAGENOME_API_KEY não está definida como variável de ambiente")

alphagenome_model = alphagenome_client.create(api_key=api_key)

alphagenome_output_metadata = alphagenome_model.output_metadata(organism=alphagenome_client.Organism.HOMO_SAPIENS)

all_dfs = []
for attr in ["rna_seq", "cage", "dnase", "atac", "chip_histone", "chip_tf",
             "splice_sites", "splice_junctions", "splice_site_usage", "procap", "contact_maps"]:
    if hasattr(alphagenome_output_metadata, attr):
        df = getattr(alphagenome_output_metadata, attr)
        df["output_type"] = attr.upper()
        all_dfs.append(df)
df_all = pd.concat(all_dfs, ignore_index=True)
df_all.to_csv("alpha_genome_all_tracks.csv", index=False)

# print(alphagenome_output_metadata.atac)  # DataFrame com todas as tracks de ATAC_SEQ

print(alphagenome_output_metadata.rna_seq)  # DataFrame com todas as tracks de RNA_SEQ
