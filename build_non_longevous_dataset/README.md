# Non-Longevous Dataset Builder

> **📁 Localização**: Este módulo está em `build_non_longevous_dataset/`

Pipeline para construir datasets de indivíduos não longevos do projeto 1000 Genomes.

## 📋 Descrição

Este programa analisa um arquivo CSV com metadados de indivíduos do projeto 1000 Genomes, permite selecionar amostras baseado em critérios personalizados, e executa análises genômicas usando `build_window_and_predict.py` para cada indivíduo selecionado.

## 🔧 Requisitos

- Python 3.8+
- Pacotes Python:
  - pandas
  - pyyaml
  - numpy
  - alphagenome (para predições)
- Ferramentas:
  - samtools
  - bcftools
- Arquivos:
  - `build_window_and_predict.py` (no diretório pai do projeto)
  - Genoma de referência GRCh38 (.fa + .fai)
  - VCFs do 1000 Genomes (filtrados e faseados)

## 📊 Formato do CSV

O arquivo CSV deve conter as seguintes colunas:

```
FamilyID,SampleID,FatherID,MotherID,Sex,Population,Superpopulation
```

Onde:
- **SampleID**: Identificador único do indivíduo (e.g., HG00096)
- **Sex**: 1 = Masculino, 2 = Feminino
- **Population**: População (e.g., ACB, GBR, CHB)
- **Superpopulation**: Superpopulação (AFR, EUR, EAS, SAS, AMR)

Exemplo:
```csv
BB01,HG01879,0,0,1,ACB,AFR
BB01,HG01880,0,0,2,ACB,AFR
Y001,HG00096,0,0,1,GBR,EUR
```

## 🚀 Uso Básico

### Passo 1: Analisar Metadados

Primeiro, analise o arquivo CSV para ver estatísticas sobre as amostras disponíveis:

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml

# OU, da raiz do projeto:
python3 build_non_longevous_dataset/build_non_longevous_dataset.py --config build_non_longevous_dataset/configs/default.yaml
```

Isso irá imprimir:
- Número total de amostras
- Quantas superpopulações existem
- Quantas pessoas em cada superpopulação
- Quantas populações em cada superpopulação
- Distribuição de sexo em cada população

### Passo 2: Configurar Seleção de Amostras

Edite o arquivo `configs/default.yaml` para configurar:

1. **Caminho do CSV**:
```yaml
data_sources:
  metadata_csv: "../../doc/1000_genomes_metadata.csv"  # Relativo ao diretório configs/
```

2. **Critérios de seleção**:
```yaml
sample_selection:
  level: "superpopulation"  # ou "population"
  samples_per_group: 2       # quantas amostras por grupo
  sex_filter: "all"          # "all", "male", ou "female"
```

3. **Parâmetros do gene a analisar**:
```yaml
build_window_params:
  gene: "CYP2B6"            # gene de interesse
  window_size: 1000000      # janela de 1 Mb
  predict: true             # executar predições AlphaGenome
  outputs: "RNA_SEQ,ATAC"   # tipos de output
  ontology: "UBERON:0002107,UBERON:0000955"  # tecidos (fígado, cérebro)
```

4. **Habilitar passos adicionais**:
```yaml
pipeline:
  steps:
    analyze_metadata: true    # Passo 1: analisar CSV
    select_samples: true      # Passo 2: selecionar amostras
    validate_vcfs: false      # Passo 3: validar VCFs (opcional)
    run_predictions: true     # Passo 4: executar predições
    generate_report: true     # Passo 5: gerar relatório
```

### Passo 3: Executar Pipeline Completo

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

## 📁 Estrutura de Saída

```
non_longevous_results/
├── metadata_statistics.json        # Estatísticas do CSV
├── selected_samples.csv            # Amostras selecionadas
├── non_longevous_dataset_checkpoint.json  # Checkpoint (idempotência)
├── processing_summary.txt          # Relatório final
├── HG00096__CYP2B6/                # Resultados por amostra/gene
│   ├── ref.window.fa
│   ├── HG00096.H1.window.fixed.fa
│   ├── HG00096.H2.window.fixed.fa
│   └── predictions_H1/
│       ├── rna_seq.npz
│       └── rna_seq_metadata.json
└── ...
```

## 🔄 Idempotência

O programa é idempotente e mantém um arquivo de checkpoint. Se a execução for interrompida:

1. Amostras já processadas **não** serão reprocessadas
2. O pipeline continuará de onde parou
3. Para reprocessar tudo, delete o arquivo de checkpoint:
   ```bash
   rm non_longevous_results/non_longevous_dataset_checkpoint.json
   ```

## ⚙️ Opções Avançadas

### Selecionar Apenas Algumas Populações

```yaml
sample_selection:
  level: "population"
  samples_per_group: 5
  include_groups: ["GBR", "CHB", "YRI"]  # apenas estas populações
```

### Excluir Populações

```yaml
sample_selection:
  exclude_groups: ["ACB", "ASW"]  # excluir estas
```

### Filtrar por Sexo

```yaml
sample_selection:
  sex_filter: "male"  # apenas masculino
```

### Desabilitar Predições AlphaGenome (mais rápido)

```yaml
build_window_params:
  predict: false  # apenas extrair sequências
```

## 📊 Exemplo de Saída (Passo 1)

```
================================================================================
ESTATÍSTICAS DO DATASET - 1000 GENOMES PROJECT
================================================================================

📊 TOTAL DE AMOSTRAS: 56

🌍 SUPERPOPULAÇÕES: 5
--------------------------------------------------------------------------------

  AFR:
    • Total de indivíduos: 16
    • Masculino: 8
    • Feminino: 8
    • Número de populações: 2
    • Populações: ACB, ASW

  AMR:
    • Total de indivíduos: 10
    • Masculino: 5
    • Feminino: 5
    • Número de populações: 2
    • Populações: MXL, PUR

  EAS:
    • Total de indivíduos: 10
    • Masculino: 5
    • Feminino: 5
    • Número de populações: 2
    • Populações: CHB, CHS

  EUR:
    • Total de indivíduos: 10
    • Masculino: 5
    • Feminino: 5
    • Número de populações: 2
    • Populações: GBR, TSI

  SAS:
    • Total de indivíduos: 10
    • Masculino: 5
    • Feminino: 5
    • Número de populações: 2
    • Populações: GIH, ITU

🏘️  POPULAÇÕES: 10
--------------------------------------------------------------------------------

  AFR:
    ACB: 10 indivíduos (♂ 5, ♀ 5)
    ASW: 6 indivíduos (♂ 3, ♀ 3)

  AMR:
    MXL: 4 indivíduos (♂ 2, ♀ 2)
    PUR: 6 indivíduos (♂ 3, ♀ 3)

  ...
```

## 🧬 Superpopulações do 1000 Genomes

- **AFR**: African (Africana)
- **AMR**: Ad Mixed American (Américas Mistas)
- **EAS**: East Asian (Leste Asiático)
- **EUR**: European (Europeia)
- **SAS**: South Asian (Sul Asiático)

## 💡 Dicas

1. **Comece com análise**: Execute apenas o passo `analyze_metadata` primeiro para entender seus dados
2. **Teste com poucos**: Use `samples_per_group: 1` ou `2` para testes rápidos
3. **Use checkpoint**: O sistema salva progresso automaticamente
4. **Ontologias múltiplas**: Separe por vírgula: `"UBERON:0002107,UBERON:0000955,CL:0002601"`
5. **VCF por cromossomo**: Certifique-se que o VCF contém o cromossomo do seu gene

## 🔍 Troubleshooting

### Erro: CSV não encontrado
```
[ERROR] Arquivo CSV não encontrado: ...
```
**Solução**: Verifique o caminho em `data_sources.metadata_csv` no YAML

### Erro: VCF pattern contém {chrom}
```
[WARN] VCF pattern contém {chrom}, mas cromossomo não foi determinado.
```
**Solução**: Forneça o caminho completo do VCF ou especifique o cromossomo

### Erro: API key não encontrado
```
RuntimeError: AlphaGenome API key not provided
```
**Solução**: 
```bash
export ALPHAGENOME_API_KEY="sua_chave_aqui"
```

## 📝 Autor

Alberto F. De Souza
Última atualização: 2025-11-04

