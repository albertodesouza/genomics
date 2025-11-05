# AlphaGenome Predictions - Guia de Uso

## 📋 Visão Geral

O script `build_window_and_predict.py` agora salva as predições completas do AlphaGenome como arrays NumPy, permitindo análises detalhadas dos dados de predição (ATAC-seq, RNA-seq, etc.) para cada nucleotídeo da sequência.

## 🚀 Executando Predições

### Exemplo básico:

```bash
# Do diretório build_non_longevous_dataset
python3 build_window_and_predict.py \
  --sample HG00096 \
  --gene CYP2B6 \
  --ref-fasta ../refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --vcf ../longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --outdir ./alphagenome_output \
  --predict \
  --outputs "ATAC" \
  --ontology "UBERON:0002107"
```

### Ver outputs disponíveis:

```bash
python3 build_window_and_predict.py --list-outputs
```

Saída exemplo:
```
Available OutputType attributes in AlphaGenome:
  ATAC
  CAGE
  DNASE
  H3K27AC
  H3K27ME3
  H3K36ME3
  H3K4ME1
  H3K4ME3
  H3K9ME3
  RNA
  ...
```

## 📁 Estrutura de Saída

Após a execução, os arquivos são organizados assim:

```
alphagenome/HG00096__CYP2B6/
├── gtf_cache.feather                         # Cache do GTF (reutilizado)
├── ref.window.fa                             # Sequência de referência (1 Mb)
├── HG00096.window.vcf.gz                     # Variantes do sample na região
├── HG00096.window.consensus_ready.vcf.gz     # Variantes filtradas
├── HG00096.H1.window.raw.fa                  # Consenso H1 (antes do ajuste)
├── HG00096.H1.window.fixed.fa                # Consenso H1 (exatos 1 milhão de bases)
├── HG00096.H2.window.raw.fa                  # Consenso H2 (antes do ajuste)
├── HG00096.H2.window.fixed.fa                # Consenso H2 (exatos 1 milhão de bases)
├── predictions_H1/                           # ⭐ Predições AlphaGenome para H1
│   ├── atac.npz                              #    Arrays NumPy (1M valores)
│   └── atac_metadata.json                    #    Metadados dos tracks
├── predictions_H2/                           # ⭐ Predições AlphaGenome para H2
│   ├── atac.npz
│   └── atac_metadata.json
├── prediction_H1.ok.txt                      # Marker de conclusão H1
└── prediction_H2.ok.txt                      # Marker de conclusão H2
```

## 📊 Analisando os Resultados

### Script de análise incluído:

```bash
# Análise básica de um arquivo
python3 ~/genomics/read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz

# Gerar plot de uma região
python3 ~/genomics/read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz \
  --plot --start 0 --end 10000 --output atac_plot.png

# Comparar haplótipos H1 vs H2
python3 ~/genomics/read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz \
  --compare alphagenome/HG00096__CYP2B6/predictions_H2/atac.npz
```

### Exemplo em Python:

```python
import numpy as np
import json
from pathlib import Path

# Carregar predições H1
data_h1 = np.load('alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz')

# Ver tracks disponíveis
print(f"Tracks: {data_h1.files}")  # Ex: ['track_0', 'track_1', ...]

# Acessar track específico
track_0 = data_h1['track_0']  # Array com ~1 milhão de valores

# Estatísticas básicas
print(f"Shape: {track_0.shape}")
print(f"Mean:  {track_0.mean():.6f}")
print(f"Std:   {track_0.std():.6f}")
print(f"Min:   {track_0.min():.6f}")
print(f"Max:   {track_0.max():.6f}")

# Carregar metadados
with open('alphagenome/HG00096__CYP2B6/predictions_H1/atac_metadata.json') as f:
    metadata = json.load(f)
    
print(f"Metadados: {metadata}")

# Analisar região específica (ex: primeiros 1000 nucleotídeos)
region = track_0[0:1000]
print(f"Média na região 0-1000: {region.mean():.6f}")

# Comparar H1 vs H2
data_h2 = np.load('alphagenome/HG00096__CYP2B6/predictions_H2/atac.npz')
track_h2 = data_h2['track_0']

# Diferença absoluta
diff = np.abs(track_0 - track_h2)
print(f"Diferença média entre H1 e H2: {diff.mean():.6f}")
print(f"Posições com diferença > 0.1: {(diff > 0.1).sum()}")

# Salvar resultados processados
np.save('diferenca_h1_h2.npy', diff)
```

## 🧬 CURIEs de Tecidos Comuns

Para usar com `--tissue`:

| CURIE | Tecido/Célula |
|-------|---------------|
| `UBERON:0002107` | Fígado (liver) |
| `UBERON:0000955` | Cérebro (brain) |
| `UBERON:0000948` | Coração (heart) |
| `UBERON:0002048` | Pulmão (lung) |
| `UBERON:0001264` | Pâncreas (pancreas) |
| `CL:0000182` | Hepatócito |
| `CL:0000540` | Neurônio |
| `CL:0000746` | Cardiomiócito |

Se não especificar `--tissue` ou usar um valor inválido, o AlphaGenome retorna predições para **todos** os tecidos/células disponíveis.

## 🔄 Idempotência

O script é completamente idempotente:

- ✅ Cache do GTF (reutilizado entre todos os genes)
- ✅ Sequências FASTA (puladas se já existem)
- ✅ VCFs processados (pulados se já existem)
- ✅ Predições (puladas se markers `.ok.txt` existem)

Você pode executar o mesmo comando múltiplas vezes e apenas os passos incompletos serão executados.

## ⚡ Performance

### Primeira execução (sem cache):
- Download GTF: ~10-30 segundos
- Extração de referência: ~1-2 segundos
- Subset VCF: ~2-5 segundos
- Consensus (H1+H2): ~3-5 segundos
- Predições AlphaGenome: ~30-60 segundos (depende da API)
- **Total: ~1-2 minutos**

### Execuções subsequentes (com cache):
- Carregamento GTF: ~0.5 segundos
- Pula todos os passos já feitos
- **Total: ~1 segundo** (se tudo já existe)

## 💾 Espaço em Disco

Por caso (sample + gene):

- Sequências FASTA: ~3-5 MB
- VCFs: ~0.5-2 MB (depende do número de variantes)
- Predições NPZ (comprimidas): ~8-20 MB por output type por haplótipo
- **Total estimado**: ~20-50 MB por caso

O cache do GTF (~50-100 MB) é compartilhado entre todos os casos.

## 🎯 Casos de Uso

### 1. Análise de impacto funcional
Compare predições entre haplótipos para identificar variantes com efeito funcional:

```python
diff = np.abs(h1_predictions - h2_predictions)
high_impact_positions = np.where(diff > threshold)[0]
```

### 2. Perfil epigenético de genes
Analise o perfil de cromatina (ATAC, H3K27ac, etc.) ao redor de um gene de interesse.

### 3. Efeitos específicos de tecido
Compare predições usando diferentes `--tissue` para ver se variantes têm efeitos tecido-específicos.

### 4. Análise populacional
Execute para múltiplos samples (ex: 1000 Genomes) e compare perfis entre populações.

## 📚 Referências

- [AlphaGenome Documentation](https://alphafold.com/alphagenome)
- [UBERON Ontology Browser](https://www.ebi.ac.uk/ols/ontologies/uberon)
- [Cell Ontology (CL)](https://www.ebi.ac.uk/ols/ontologies/cl)

## 🐛 Troubleshooting

### Erro: "Invalid ontology_curie"
Use CURIEs no formato `TIPO:ID` (ex: `UBERON:0002107`), não texto livre.

### Erro: "Output type not found"
Use `--list-outputs` para ver nomes válidos. Use o nome exato (ex: `ATAC`, não `ATAC-seq`).

### Arrays vazios ou None
Algumas combinações output/tissue podem não ter dados. Verifique os warnings no log.

### Falta memória
As predições usam ~2-4 GB de RAM. Para múltiplos samples, processe sequencialmente.

