# verify_processed_dataset.py - Verificação de Dataset Processado

## 📋 Visão Geral

Programa para verificar e comparar dados processados de sequências genômicas com predições do AlphaGenome. Visualiza tracks (genes × ontologias RNA-seq) em gráficos de linha sobrepostos para detectar possíveis bugs no pipeline de processamento e validar consistência dos dados.

## 🎯 Objetivo

Detectar inconsistências e validar comparações entre:
- **Cache processado**: Dados normalizados armazenados em `.pt` para treinamento
- **Dataset original**: Dados processados em `.npz` a partir de genomas individuais
- **AlphaGenome (referência)**: Predições brutas da API usando genoma de referência
- **AlphaGenome (individual)**: Predições brutas da API usando genoma do indivíduo

## ✨ Funcionalidades

### 🔄 Modos de Comparação

O programa suporta 4 modos de comparação distintos:

1. **`alphagenome_ref_x_dataset_dir`**: 
   - Compara predições do AlphaGenome (genoma de referência) vs. dados do dataset processado
   - Útil para validar se o dataset reflete corretamente o genoma de referência

2. **`alphagenome_ind_x_dataset_dir`**: 
   - Compara predições do AlphaGenome (genoma individual com variantes) vs. dados do dataset processado
   - Valida se o processamento de variantes individuais está correto
   - Usa `build_window_and_predict.py` como biblioteca

3. **`dataset_dir_x_cache_dir`**: 
   - Compara dados do dataset original (`.npz`) vs. cache processado (`.pt`)
   - Valida pipeline de normalização e transformação

4. **`alphagenome_x_alphagenome_ref`**: 
   - Compara duas formas de chamar a API do AlphaGenome
   - `predict_interval` (sem FASTA) vs. `predict_sequence` (com FASTA extraído)
   - Valida consistência da API

### 🎮 Modos de Navegação Interativa

#### Modo "single" (padrão)
Navegação tradicional por um indivíduo por vez:
- **← (seta esquerda)**: Retrocede para amostra anterior
- **→ (seta direita)**: Avança para próxima amostra
- **Q**: Sai do programa

#### Modo "comparison" (novo!)
Comparação interativa entre dois indivíduos:
- **← → (setas)**: Avança/retrocede ambos os indivíduos simultaneamente
- **A**: Retrocede apenas o segundo indivíduo
- **D**: Avança apenas o segundo indivíduo
- **W**: Avança para o próximo gene (ambos indivíduos)
- **Z**: Retrocede para o gene anterior (ambos indivíduos)
- **Q**: Sai do programa

No modo "comparison":
- Visualiza 6 tracks empilhadas verticalmente (3 ontologias × 2 strands)
- Primeiro indivíduo em azul sólido
- Segundo indivíduo em vermelho tracejado
- Legenda mostra: `sample_id (population/superpopulation)`
- Gene atual exibido no título principal

### 🧬 Filtro de Genes
- Visualizar todos os 11 genes (66 tracks)
- Visualizar apenas um gene específico (6 tracks)
- Visualizar múltiplos genes selecionados

### ⚙️ Configuração via YAML
- Todas as configurações em arquivo YAML
- Múltiplas configurações predefinidas
- Fácil customização

## 🚀 Uso

### Uso Básico (modo interativo, todos os genes)

```bash
cd neural_ancestry_predictor

python3 verify_processed_dataset.py --config configs/verify_processed_dataset.yaml
```

### Verificar apenas um gene específico (ex: MC1R)

```bash
python3 verify_processed_dataset.py --config configs/verify_tyr_only.yaml
```

### Modo de Comparação entre Dois Indivíduos

```yaml
# configs/verify_comparison.yaml
interactive_mode: true
interactive_comparison_mode: "comparison"  # Ativa modo de comparação
comparison_mode: "dataset_dir_x_cache_dir"
gene_filter: "MC1R"  # Recomendado: um gene por vez
```

```bash
python3 verify_processed_dataset.py --config configs/verify_comparison.yaml
# Use ← → para navegar ambos
# Use A D para navegar apenas o segundo indivíduo
# Use W Z para mudar o gene
```

### Criar configuração customizada

Crie um arquivo `.yaml` no diretório `configs/`:

```yaml
# configs/my_verification.yaml
cache_dir: "/path/to/cache"
dataset_dir: "/path/to/dataset"
split: "test"
index: 0
gene_filter: "TYR"  # null para todos, "GENE" para um, ["G1", "G2"] para múltiplos

# Modo de comparação
comparison_mode: "alphagenome_ref_x_dataset_dir"

# Modo de navegação
interactive_mode: true
interactive_comparison_mode: "single"  # ou "comparison"

show_navigation_help: true
verbose_metrics: true
show_stats_in_plot: true
save_plots: false
output_dir: null
```

Depois execute:

```bash
python3 verify_processed_dataset.py --config configs/my_verification.yaml
```

## 📝 Configurações Disponíveis

### Arquivo de Configuração YAML

```yaml
# ═══════════════════════════════════════════════════════════════
# DIRETÓRIOS (obrigatórios)
# ═══════════════════════════════════════════════════════════════

cache_dir: "/path/to/cache"           # Cache processado (.pt, metadata.json)
dataset_dir: "/path/to/dataset"       # Dataset original (.npz)

# ═══════════════════════════════════════════════════════════════
# SELEÇÃO DE AMOSTRA
# ═══════════════════════════════════════════════════════════════

split: "test"                         # train, val ou test
index: 0                              # Índice inicial (muda com ← →)
sample_id: null                       # ID específico (opcional, sobrescreve split/index)

# ═══════════════════════════════════════════════════════════════
# MODO DE COMPARAÇÃO
# ═══════════════════════════════════════════════════════════════

comparison_mode: "dataset_dir_x_cache_dir"
# Opções:
#   - alphagenome_ref_x_dataset_dir: AlphaGenome (ref) vs dataset
#   - alphagenome_ind_x_dataset_dir: AlphaGenome (ind) vs dataset
#   - dataset_dir_x_cache_dir: Dataset vs cache (padrão)
#   - alphagenome_x_alphagenome_ref: API interval vs sequence

# ═══════════════════════════════════════════════════════════════
# FILTRO DE GENES
# ═══════════════════════════════════════════════════════════════

gene_filter: null                     # Opções:
                                      # null: todos os 11 genes (66 tracks)
                                      # "TYR": apenas gene TYR (6 tracks)
                                      # ["TYR", "TYRP1"]: múltiplos genes (12 tracks)

# ═══════════════════════════════════════════════════════════════
# NAVEGAÇÃO INTERATIVA
# ═══════════════════════════════════════════════════════════════

interactive_mode: true                # true: navegação com teclado
                                      # false: visualiza apenas uma amostra

interactive_comparison_mode: "single" # "single": um indivíduo por vez (← →)
                                      # "comparison": dois indivíduos (← → A D W Z)

show_navigation_help: true            # Mostra instruções no gráfico

# ═══════════════════════════════════════════════════════════════
# ALPHAGENOME API (OPCIONAL)
# ═══════════════════════════════════════════════════════════════

alphagenome_api:
  enabled: false                      # Habilita chamadas à API
  api_key: null                       # Usa ALPHAGENOME_API_KEY do ambiente
  rate_limit_delay: 0.5               # Delay entre chamadas (segundos)
  ontology_terms: ["CL:1000458", "CL:0000346", "CL:2000092"]

# ═══════════════════════════════════════════════════════════════
# RAW MODE (OPCIONAL)
# ═══════════════════════════════════════════════════════════════

raw_mode:
  enabled: false                      # Modo raw (apenas AlphaGenome, sem cache)
  source: "files"                     # "files" (.npz) ou "api" (chamada API)
  window_size_key: "SEQUENCE_LENGTH_16KB"  
  # Opções: SEQUENCE_LENGTH_2KB, SEQUENCE_LENGTH_16KB, 
  #         SEQUENCE_LENGTH_100KB (128 KiB), SEQUENCE_LENGTH_500KB (512 KiB),
  #         SEQUENCE_LENGTH_1MB

# ═══════════════════════════════════════════════════════════════
# VISUALIZAÇÃO
# ═══════════════════════════════════════════════════════════════

verbose_metrics: true                 # Exibir métricas detalhadas no console
show_stats_in_plot: true              # Mostrar estatísticas no gráfico

# ═══════════════════════════════════════════════════════════════
# SALVAMENTO (opcional)
# ═══════════════════════════════════════════════════════════════

save_plots: false                     # Salvar gráficos automaticamente
output_dir: null                      # Diretório para salvar (cria se não existe)
output_prefix: "verify"               # Prefixo dos arquivos (ex: verify_HG00120.png)
```

## 🌐 Modo API do AlphaGenome

### O que é?

O modo API permite verificar o dataset comparando os dados do cache diretamente com predições **em tempo real** da API do AlphaGenome, ao invés de usar os arquivos `.npz` pré-computados.

### Quando usar?

- ✅ Verificar se os `.npz` originais foram gerados corretamente
- ✅ Testar com sequências personalizadas (variantes individuais)
- ✅ Validação end-to-end completa do pipeline
- ✅ Debugging de inconsistências nos dados processados

### Como ativar?

1. Configure a API key:
```bash
export ALPHAGENOME_API_KEY="your_api_key_here"
```

2. Crie/edite configuração YAML:
```yaml
# Habilitar modo API
alphagenome_api:
  enabled: true  # Ativa chamadas à API
  api_key: null  # Usa ALPHAGENOME_API_KEY do ambiente
  rate_limit_delay: 0.5  # Delay entre chamadas (segundos)
  ontology_terms: ["CL:1000458", "CL:0000346", "CL:2000092"]

# Recomendado: testar um gene por vez
gene_filter: "MC1R"
```

3. Execute:
```bash
python3 verify_processed_dataset.py --config configs/verify_api_test.yaml
```

### Requisitos

- 📦 Pacote `alphagenome` instalado (`pip install alphagenome`)
- 🔑 API key válida do AlphaGenome
- 🌐 Conexão com internet
- 📁 Arquivos `.fa` das sequências em `dataset_dir/individuals/{sample}/windows/{gene}/`

### Importante

⚠️ **Quota de API**: Cada gene consome 1 chamada à API. Use `gene_filter` para economizar!  
⚠️ **Rate limiting**: Configure `rate_limit_delay` adequadamente (padrão: 0.5s)  
⚠️ **Ontologias**: Devem corresponder **exatamente** às usadas na criação do dataset original

### Constantes AlphaGenome

O AlphaGenome usa constantes específicas para tamanhos de janela:

| Constante | Tamanho (bp) | Tamanho (KiB) | Uso |
|-----------|--------------|---------------|-----|
| `SEQUENCE_LENGTH_2KB` | 2048 | 2 KiB | Testes rápidos |
| `SEQUENCE_LENGTH_16KB` | 16384 | 16 KiB | Genes pequenos |
| `SEQUENCE_LENGTH_100KB` | 131072 | 128 KiB | Genes médios |
| `SEQUENCE_LENGTH_500KB` | 524288 | 512 KiB | Genes grandes |
| `SEQUENCE_LENGTH_1MB` | 1048576 | 1 MiB | Regiões extensas |

⚠️ **Nota**: Os nomes das constantes usam KB (base 1000), mas os tamanhos são em potências de 2 (KiB).

## 📊 Saída

### Console (Modo "single")

```
════════════════════════════════════════════════════════════
       VERIFICAÇÃO DE DATASET PROCESSADO                   
════════════════════════════════════════════════════════════

Modo interativo ativado:
  • Split: test
  • Total de amostras: 13
  • Índice inicial: 0
  • Use ← → para navegar, 'q' para sair

✓ Sample: NA19472 (índice 0, global 65)

═══════════════════════════════════════════════════════
           MÉTRICAS DE COMPARAÇÃO                      
═══════════════════════════════════════════════════════
┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┓
┃ Métrica               ┃              Valor ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━┩
│ MAE Média (global)    │           0.004782 │
│ MAE Máximo            │           0.048390 │
│ MAE Mínimo            │           0.000001 │
│ Track com maior MAE   │ 45 (gene 7, ont 3) │
│                       │                    │
│ Correlação Média      │           0.337512 │
│ Correlação Mínima     │          -0.043882 │
│ Correlação Máxima     │           0.988911 │
│ Track com menor corr. │ 55 (gene 9, ont 1) │
└───────────────────────┴────────────────────┘
═══════════════════════════════════════════════════════
```

### Console (Modo "comparison")

```
════════════════════════════════════════════════════════════
       VERIFICAÇÃO DE DATASET PROCESSADO                   
════════════════════════════════════════════════════════════

Modo de comparação interativa ativado:
  • Split: train
  • Total de amostras: 54
  • Total de genes: 11
  • Use ← → (ambos), A D (ind2), W Z (genes), Q (sair)

✓ Mostrando: MC1R
  Individual 1: HG00120 (GBR/EUR)
  Individual 2: HG00120 (GBR/EUR)
```

### Gráfico (Modo "single")

O gráfico mostra:
- **Azul sólido**: Dados do primeiro conjunto (cache, dataset, ou AlphaGenome)
- **Vermelho tracejado**: Dados do segundo conjunto (AlphaGenome, cache, ou referência)
- **Linhas cinza pontilhadas**: Separação entre genes
- **Labels do eixo Y**: Ontologia e tipo celular em duas linhas
  - Linha 1: `CL:XXXXX (±)` (código da ontologia e strand)
  - Linha 2: Nome do tipo celular (ex: `Melanocyte`)
- **Box com estatísticas**: MAE e correlação por track
- **Instruções de navegação**: Canto inferior central (se modo interativo)

### Gráfico (Modo "comparison")

O gráfico mostra:
- **6 subplots empilhados verticalmente** (3 ontologias × 2 strands)
- **Azul sólido**: Tracks do primeiro indivíduo
- **Vermelho tracejado**: Tracks do segundo indivíduo
- **Título principal**: Gene atual e informação dos dois indivíduos
- **Labels do eixo Y**: 
  - `CL:1000458 (+)` / `Melanocyte`
  - `CL:0000346 (+)` / `Dermal Papilla`
  - `CL:2000092 (+)` / `Keratinocyte`
  - (e versões (-) para strand negativo)
- **Legenda**: `sample_id (population/superpopulation)` para cada indivíduo

### Ontologias RNA-seq

O dataset usa 3 ontologias de tipos celulares, cada uma com 2 strands (+/-):

| Código CL | Tipo Celular | Descrição |
|-----------|--------------|-----------|
| CL:1000458 | Melanocyte | Células produtoras de melanina |
| CL:0000346 | Dermal Papilla | Células da papila dérmica |
| CL:2000092 | Keratinocyte | Células da epiderme |

Cada gene tem **6 tracks**: 3 ontologias × 2 strands (+ e -)

## 🧬 Genes Disponíveis

1. **SLC24A5** - Solute carrier family 24 member 5
2. **SLC45A2** - Solute carrier family 45 member 2
3. **OCA2** - OCA2 melanosomal transmembrane protein
4. **HERC2** - HECT and RLD domain containing E3 ubiquitin protein ligase 2
5. **MC1R** - Melanocortin 1 receptor
6. **EDAR** - Ectodysplasin A receptor
7. **MFSD12** - Major facilitator superfamily domain containing 12
8. **DDB1** - Damage specific DNA binding protein 1
9. **TCHH** - Trichohyalin
10. **TYR** - Tyrosinase
11. **TYRP1** - Tyrosinase related protein 1

Cada gene tem 6 tracks correspondentes a 3 ontologias × 2 strands.

## 🔍 Interpretação dos Resultados

### MAE (Mean Absolute Error)
- **< 0.01**: ✅ Excelente correspondência
- **0.01 - 0.05**: ⚠️ Boa correspondência, pequenas diferenças
- **> 0.05**: ❌ Possível bug no pipeline

### Correlação de Pearson
- **> 0.9**: ✅ Excelente correlação
- **0.7 - 0.9**: ✅ Boa correlação
- **0.5 - 0.7**: ⚠️ Correlação moderada
- **< 0.5**: ❌ Correlação fraca - investigar

### No Gráfico
1. **Sobreposição**: Linhas azul e vermelha devem estar próximas
2. **Padrões**: Formas das curvas devem ser similares
3. **Outliers**: Tracks com grande divergência precisam de investigação

### Modo "comparison"
- No modo de comparação, não há métricas MAE/correlação
- O foco é na **comparação visual** das tracks de dois indivíduos
- Diferenças entre indivíduos podem indicar variações genéticas reais

## 📚 Exemplos de Uso

### Exemplo 1: Percorrer todas as amostras do teste

```yaml
# configs/scan_test_all.yaml
cache_dir: "/dados/GENOMICS_DATA/top3/non_longevous_results_runs_genes/datasets/rna_seq_H1_100000_ds1_log_split0.7-0.15-0.15_seed5"
dataset_dir: "/dados/GENOMICS_DATA/top3/non_longevous_results_genes"
split: "test"
index: 0
gene_filter: null  # Todos os genes
interactive_mode: true
interactive_comparison_mode: "single"
comparison_mode: "dataset_dir_x_cache_dir"
```

```bash
python3 verify_processed_dataset.py --config configs/scan_test_all.yaml
# Pressione → para navegar pelas 13 amostras do teste
# Pressione 'q' para sair
```

### Exemplo 2: Verificar gene específico em todas as amostras

```yaml
# configs/verify_tyr_all.yaml
gene_filter: "TYR"
interactive_mode: true
interactive_comparison_mode: "single"
split: "train"  # 54 amostras
comparison_mode: "dataset_dir_x_cache_dir"
```

```bash
python3 verify_processed_dataset.py --config configs/verify_tyr_all.yaml
# Navegue com → para ver TYR em todas as amostras de treino
```

### Exemplo 3: Comparar dois indivíduos interativamente

```yaml
# configs/compare_individuals.yaml
interactive_mode: true
interactive_comparison_mode: "comparison"
comparison_mode: "dataset_dir_x_cache_dir"
gene_filter: "MC1R"  # Recomendado: um gene por vez
split: "train"
```

```bash
python3 verify_processed_dataset.py --config configs/compare_individuals.yaml
# ← → : navega ambos os indivíduos
# A D : navega apenas o segundo indivíduo
# W Z : muda o gene sendo exibido
# Q   : sai
```

### Exemplo 4: Validar AlphaGenome vs Dataset (referência)

```yaml
# configs/validate_alphagenome_ref.yaml
comparison_mode: "alphagenome_ref_x_dataset_dir"
gene_filter: "MC1R"
interactive_mode: true
interactive_comparison_mode: "single"
sample_id: "HG02445"
```

```bash
python3 verify_processed_dataset.py --config configs/validate_alphagenome_ref.yaml
# Compara AlphaGenome (genoma de referência) com dataset processado
```

### Exemplo 5: Validar AlphaGenome vs Dataset (individual)

```yaml
# configs/validate_alphagenome_ind.yaml
comparison_mode: "alphagenome_ind_x_dataset_dir"
gene_filter: "MC1R"
interactive_mode: false
sample_id: "HG02445"
```

```bash
python3 verify_processed_dataset.py --config configs/validate_alphagenome_ind.yaml
# Compara AlphaGenome (genoma individual) com dataset processado
# Usa build_window_and_predict.py para gerar predições individuais
```

### Exemplo 6: Salvar gráficos de todas as amostras

```yaml
# configs/save_all_plots.yaml
interactive_mode: false
save_plots: true
output_dir: "verification_plots"
output_prefix: "verify"
gene_filter: "TYR"
comparison_mode: "dataset_dir_x_cache_dir"
```

```bash
# Criar script para processar todas as amostras
for i in {0..12}; do
    sed "s/index: 0/index: $i/" configs/save_all_plots.yaml > /tmp/config_$i.yaml
    python3 verify_processed_dataset.py --config /tmp/config_$i.yaml
done
```

## 🛠️ Workflow Recomendado

### 1. Verificação Inicial (Visão Geral)
```bash
# Verificar algumas amostras com todos os genes
python3 verify_processed_dataset.py --config configs/verify_processed_dataset.yaml
# Use → para ver 3-5 amostras diferentes
# Se MAE > 0.05 ou Corr < 0.5 → investigar
```

### 2. Investigação por Gene
```bash
# Se problema detectado, isolar gene problemático
python3 verify_processed_dataset.py --config configs/verify_tyr_only.yaml
# Verificar se problema é específico do gene ou geral
```

### 3. Validação AlphaGenome
```bash
# Validar dados contra AlphaGenome (referência)
# Editar config para comparison_mode: "alphagenome_ref_x_dataset_dir"
python3 verify_processed_dataset.py --config configs/verify_alphagenome_ref.yaml

# Validar dados contra AlphaGenome (individual)
# Editar config para comparison_mode: "alphagenome_ind_x_dataset_dir"
python3 verify_processed_dataset.py --config configs/verify_alphagenome_ind.yaml
```

### 4. Comparação entre Indivíduos
```bash
# Comparar padrões de expressão entre indivíduos
# Editar config para interactive_comparison_mode: "comparison"
python3 verify_processed_dataset.py --config configs/compare_individuals.yaml
# Use A D para navegar o segundo indivíduo independentemente
# Use W Z para mudar genes
```

### 5. Varredura Sistemática
```bash
# Verificar todos os splits
python3 verify_processed_dataset.py --config configs/verify_train.yaml
python3 verify_processed_dataset.py --config configs/verify_val.yaml  
python3 verify_processed_dataset.py --config configs/verify_test.yaml
```

### 6. Documentação
```bash
# Salvar gráficos para relatório
# Configure save_plots: true e output_dir
python3 verify_processed_dataset.py --config configs/save_for_report.yaml
```

## 🐛 Troubleshooting

### Erro: "Arquivo de configuração não existe"
- Verificar caminho do arquivo YAML
- Verificar que está no diretório correto

### Erro: "Cache dir não existe"
- Verificar `cache_dir` no YAML
- Verificar permissões de leitura

### Erro: "Gene não encontrado"
- Verificar ortografia do nome do gene
- Usar lista completa: SLC24A5, SLC45A2, OCA2, HERC2, MC1R, EDAR, MFSD12, DDB1, TCHH, TYR, TYRP1

### Erro: "Invalid window_size_key: SEQUENCE_LENGTH_512KB"
- AlphaGenome usa `SEQUENCE_LENGTH_500KB` para 524288 bp (512 KiB)
- AlphaGenome usa `SEQUENCE_LENGTH_100KB` para 131072 bp (128 KiB)
- Verificar constantes corretas na tabela acima

### Erro: "Sequence length X not supported by the model"
- AlphaGenome suporta apenas: 2048, 16384, 131072, 524288, 1048576 bp
- Verificar se o dataset foi gerado com um desses tamanhos
- Usar `window_size_key` correspondente ao tamanho do dataset

### MAE muito alto (> 0.05)
1. Verificar que cache foi gerado com mesmos parâmetros
2. Verificar método de normalização
3. Verificar se dados do AlphaGenome estão corretos
4. Usar modo `alphagenome_ref_x_dataset_dir` para validar

### Janela não responde a teclas
- Garantir que janela matplotlib está em foco
- Em alguns sistemas, clicar na janela antes de pressionar teclas

### Gráfico muito "cheio"
- Usar `gene_filter` para visualizar menos tracks
- Exemplo: `gene_filter: "TYR"` mostra apenas 6 tracks
- No modo "comparison", usar um gene por vez

### Modo alphagenome_ind_x_dataset_dir lento
- Este modo chama `build_window_and_predict.py` para gerar predições
- Usa arquivos temporários em `/tmp/GENOMICS_DATA/top3`
- Requer acesso ao FASTA de referência e VCF
- Limpa arquivos temporários automaticamente ao final

## 🔧 Detalhes Técnicos

### Window Size e Context

O programa implementa lógica sofisticada para lidar com tamanhos de janela:

1. **Detecção automática**: Lê o dataset e detecta seu tamanho original (full_length)
2. **Predição AlphaGenome**: Chama API com o mesmo full_length para manter contexto
3. **Visualização**: Extrai janela central menor (viz_length) de ambos para comparação
4. **Centralização**: Usa função `extract_center_window()` consistente

Exemplo:
- Dataset gerado com 524288 bp (512 KiB)
- AlphaGenome chamado com 524288 bp (mesmo contexto)
- Visualização mostra apenas 16384 bp (centro)
- Ambos extraídos do mesmo centro → alinhamento perfeito

### Boundary Conditions

Genes próximos ao início/fim de cromossomos são tratados corretamente:

1. **Coordenadas clipped**: `samtools faidx` pode retornar sequências menores
2. **Padding com 'N'**: Sequências são padded no início/fim conforme necessário
3. **Alinhamento mantido**: Posições relativas preservadas mesmo com padding

### Centralização de Lógica

Função `extract_center_window()` usada consistentemente:
- Remove duplicação de código
- Garante mesma lógica de centralização
- Corrige bug de off-by-one em versões anteriores

## 📦 Dependências

- Python 3.10+
- PyTorch
- NumPy
- Matplotlib
- SciPy
- PyYAML
- Rich
- Pandas
- alphagenome (opcional, para modo API)
- samtools (opcional, para modo alphagenome_ind)
- bcftools (opcional, para modo alphagenome_ind)

Todas já instaladas no ambiente `genomics`.

## 📄 Arquivos de Exemplo

O repositório inclui:
- `configs/verify_processed_dataset.yaml` - Configuração padrão (todos os genes, modo single)
- `configs/verify_tyr_only.yaml` - Apenas gene MC1R (exemplo de filtro)
- `configs/verify_raw_test.yaml` - Modo raw (AlphaGenome apenas)

## 🔄 Histórico de Versões

| Data | Versão | Mudanças |
|------|--------|----------|
| 2025-11-23 | 1.0 | Versão inicial com navegação interativa |
| 2025-11-24 | 1.1 | Adicionado filtro por gene e modo API |
| 2025-11-25 | 2.0 | Adicionados comparison modes e interactive_comparison_mode |
| 2025-11-25 | 2.1 | Correção de window_size, boundary conditions, AlphaGenome constants |
| 2025-11-25 | 2.2 | Labels de ontologia em duas linhas, modo "comparison" refinado |

## 👥 Autor

ChatGPT (for Alberto)  
Created: 2025-11-23  
Updated: 2025-11-25

## 📝 Notas Finais

Este programa é uma ferramenta essencial para garantir a qualidade e consistência do pipeline de processamento de dados genômicos. Use-o regularmente durante o desenvolvimento e antes de treinar modelos para evitar treinar com dados incorretos.

Para questões técnicas ou bugs, consulte a documentação do código ou entre em contato.
