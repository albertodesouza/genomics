# 🚀 Início Rápido: Non-Longevous Dataset Builder

## ⚡ Teste Rápido (5 minutos)

### Passo 1: Analisar Metadados

Execute o programa com o CSV configurado:

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

Você verá estatísticas sobre:
- 5 superpopulações (AFR, AMR, EAS, EUR, SAS)
- 10 populações
- 56 indivíduos totais
- Distribuição por sexo em cada população

### Passo 2: Configurar Para Seu Projeto

1. **Prepare seu CSV com metadados do 1000 Genomes**:
   - Baixe de: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/
   - Ou use o exemplo fornecido para testes

2. **Edite `configs/default.yaml`**:

```yaml
data_sources:
  # Seu CSV com metadados (caminho relativo ao diretório configs/)
  metadata_csv: "../../doc/seu_arquivo.csv"
  
  # Referência GRCh38 (caminho relativo ao diretório configs/)
  reference:
    fasta: "../../caminho/para/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  
  # Padrão de VCFs (por cromossomo)
  vcf_pattern: "/caminho/para/vcfs/1kGP_high_coverage.{chrom}.vcf.gz"

sample_selection:
  # Escolha: "superpopulation" ou "population"
  level: "superpopulation"
  
  # Quantas amostras por grupo
  samples_per_group: 2
  
  # Filtro de sexo: "all", "male", ou "female"
  sex_filter: "all"

build_window_params:
  # Gene de interesse
  gene: "CYP2B6"
  
  # Executar predições AlphaGenome
  predict: true
  
  # Tipos de output
  outputs: "RNA_SEQ,ATAC"
  
  # Tecidos/células (CURIEs)
  ontology: "UBERON:0002107,UBERON:0000955"  # fígado, cérebro

pipeline:
  steps:
    analyze_metadata: true    # ✓ Análise
    select_samples: true      # ✓ Seleção
    run_predictions: true     # ✓ Executar build_window_and_predict.py
    generate_report: true     # ✓ Relatório final
```

3. **Configure a API Key do AlphaGenome** (se usar predições):

```bash
export ALPHAGENOME_API_KEY="sua_chave_aqui"
```

### Passo 3: Executar Pipeline Completo

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

## 📊 Saídas Esperadas

```
non_longevous_results/
├── metadata_statistics.json              # Estatísticas do CSV
├── selected_samples.csv                  # Amostras selecionadas
├── non_longevous_dataset_checkpoint.json # Checkpoint (idempotência)
├── processing_summary.txt                # Relatório final
│
├── HG00096__CYP2B6/                      # Resultados por amostra
│   ├── ref.window.fa                     # Referência
│   ├── HG00096.H1.window.fixed.fa        # Haplótipo 1
│   ├── HG00096.H2.window.fixed.fa        # Haplótipo 2
│   ├── predictions_H1/                   # Predições H1
│   │   ├── rna_seq.npz
│   │   ├── rna_seq_metadata.json
│   │   ├── atac.npz
│   │   └── atac_metadata.json
│   └── predictions_H2/                   # Predições H2
│       └── ...
│
└── HG00097__CYP2B6/                      # Próxima amostra
    └── ...
```

## 🎯 Casos de Uso Comuns

### Caso 1: Comparar Populações Africana e Europeia

```yaml
sample_selection:
  level: "superpopulation"
  samples_per_group: 10
  include_groups: ["AFR", "EUR"]
```

### Caso 2: Apenas Mulheres de Populações Específicas

```yaml
sample_selection:
  level: "population"
  samples_per_group: 5
  include_groups: ["GBR", "CHB", "YRI"]
  sex_filter: "female"
```

### Caso 3: Todas as Superpopulações (Balanceado)

```yaml
sample_selection:
  level: "superpopulation"
  samples_per_group: 20
  include_groups: []  # todas
  sex_filter: "all"
```

### Caso 4: Apenas Extrair Sequências (Sem Predições)

```yaml
build_window_params:
  gene: "BRCA1"
  predict: false       # Desabilita AlphaGenome
  skip_h2: false       # Mantém ambos os haplótipos
```

## 🔧 Resolução de Problemas

### Erro: CSV não encontrado

```
[ERROR] Arquivo CSV não encontrado
```

**Solução**: Verifique o caminho em `data_sources.metadata_csv`

### Erro: VCF pattern com {chrom}

```
[WARN] VCF pattern contém {chrom}, mas cromossomo não foi determinado
```

**Solução**: Atualmente, o programa assume que você conhece qual cromossomo contém seu gene. Para o gene CYP2B6 (cromossomo 19), configure:

```yaml
vcf_pattern: "/caminho/para/1kGP_high_coverage.chr19.vcf.gz"
```

Alternativamente, deixe o `build_window_and_predict.py` determinar automaticamente baixando o GTF.

### Erro: AlphaGenome API Key

```
RuntimeError: AlphaGenome API key not provided
```

**Solução**:
```bash
export ALPHAGENOME_API_KEY="sua_chave"
# Ou adicione ao ~/.bashrc
```

## 💡 Dicas de Performance

1. **Comece pequeno**: Use `samples_per_group: 1` ou `2` para testes
2. **Desabilite H2**: Use `skip_h2: true` para processar 2x mais rápido
3. **Menos ontologias**: Use 1-3 tecidos específicos em vez de todos
4. **Paralelização**: Configure `n_workers: 8` (ajuste para seu CPU)
5. **Checkpoint**: Se interrompido, o programa continua de onde parou

## 📚 Documentação Completa

Veja `README_NON_LONGEVOUS_DATASET.md` para documentação detalhada.

## 🧬 Genes Comuns para Análise

- **CYP2B6** (chr19): Metabolismo de drogas
- **APOE** (chr19): Risco de Alzheimer
- **BRCA1** (chr17): Câncer de mama
- **BRCA2** (chr13): Câncer de mama/ovário
- **FOXO3** (chr6): Longevidade
- **TP53** (chr17): Supressor tumoral

## 🎓 Exemplo Completo de Workflow

```bash
# 1. Entrar no diretório do módulo
cd build_non_longevous_dataset

# 2. Analisar dados disponíveis
python3 build_non_longevous_dataset.py --config configs/default.yaml

# 3. Editar configuração baseado nas estatísticas
nano configs/default.yaml

# 4. Habilitar steps adicionais (editar YAML)
# select_samples: true
# run_predictions: true
# generate_report: true

# 5. Executar pipeline completo
python3 build_non_longevous_dataset.py --config configs/default.yaml

# 6. Verificar resultados
ls -lh ../non_longevous_results/
cat ../non_longevous_results/processing_summary.txt

# 7. Se interrompido, simplesmente execute novamente
# O checkpoint garantirá que amostras já processadas sejam puladas
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

## ✅ Checklist de Preparação

Antes de executar o pipeline completo:

- [ ] CSV com metadados preparado
- [ ] Genoma de referência GRCh38 (.fa + .fai) disponível
- [ ] VCFs do 1000 Genomes baixados
- [ ] VCFs indexados (.tbi)
- [ ] AlphaGenome API key configurada (se usar predições)
- [ ] Espaço em disco suficiente (~500MB-2GB por amostra)
- [ ] Configuração YAML revisada e ajustada

---

**Pronto!** Você está preparado para construir seu dataset de não longevos! 🎉

