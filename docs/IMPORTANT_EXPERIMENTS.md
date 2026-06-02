# Important Experiments

Este documento registra experimentos importantes de pesquisa para separar resultados consolidados dos muitos configs exploratorios de desenvolvimento.

Use este arquivo para:

- guardar resultados historicos relevantes;
- mapear cada resultado para o config/comando atual;
- registrar gaps de reprodutibilidade;
- adicionar novos resultados importantes com contexto suficiente para repeticao.

## Dataset E Convenções

- Dataset canonico para novos experimentos genotype: `/dados/GENOMICS_DATA/v1/1kG_high_coverage`.
- Runs genotype: `results/genotype_based_predictor/runs/`.
- Cache genotype: `results/cache/genotype_based_predictor/`.
- Comandos genotype devem usar `python3 -m genomics_cli genotype ...`.
- Comandos MLC/SNP ainda rodam via `snp_ancestry_predictor/snp_ancestry_predictor.py`.
- `top3` esta descontinuado para novas analises, mas alguns configs MLC/SNP historicos ainda apontam para paths legados e precisam migracao futura.

Validacoes recomendadas antes de reproduzir:

```bash
python3 -m genomics_cli audit-configs --fail-on-active-legacy
python3 -m genomics_cli audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 3 --fail-on-missing
```

## Resultados Historicos

Resultados no conjunto de teste antes da refatoracao de unificacao dos pipelines:

| Experimento | Haps | Precision | Recall | F1 | Accuracy | Config/comando atual | Status |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| MLC SNP genome-wide panel | 2 | 0.97 | 0.97 | 0.97 | 0.97 | derivar de `snp_ancestry_predictor/configs/default.yaml` com painel genome-wide e `prediction.haplotype_mode: H1+H2` | precisa config canonico |
| MLC SNP gene-window panel | 2 | 0.72 | 0.71 | 0.70 | 0.71 | derivar de `snp_ancestry_predictor/configs/default.yaml` com `panel_23andme_V5_genes_1000_w32768.txt` e `prediction.haplotype_mode: H1+H2` | precisa config canonico |
| PCA+RF AlphaGenome raw center crop | 1 | 0.54 | 0.65 | 0.58 | 0.65 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_pca400_rf.yaml` | pronto, requer sklearn/joblib |
| PCA+XGBoost AlphaGenome raw center crop | 1 | 0.73 | 0.74 | 0.73 | 0.74 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_pca400_xgboost.yaml` | pronto, requer xgboost |
| CNN AlphaGenome raw center crop | 1 | 0.75 | 0.74 | 0.74 | 0.74 | `genotype_based_predictor/configs/neural_legacy/default.yaml` ou config equivalente de 3 ontologias/H1 | aproximado |
| CNN AlphaGenome raw center crop | 2 | 0.75 | 0.76 | 0.75 | 0.76 | `genotype_based_predictor/configs/neural_legacy/genes_1000_all.yaml` | pronto para baseline legado |
| CNN DITA AG+MI+MD+MV | 2 | 0.91 | 0.91 | 0.91 | 0.91 | `genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml` | pronto |
| CNN DITA MI+MD+MV | 2 | 0.92 | 0.92 | 0.92 | 0.92 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_masks_only.yaml` | pronto |
| CNN DITA MI+MD+MV+MS | 2 | 0.91 | 0.91 | 0.91 | 0.91 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_masks_snp.yaml` | pronto |
| CNN DITA AG0+MI+MD+MV+MS | 2 | 0.91 | 0.91 | 0.91 | 0.91 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_variant_signal_mask.yaml` | pronto |
| CNN DITA AGN+MI+MD+MV+MS | 2 | 0.90 | 0.90 | 0.90 | 0.90 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_delta_reference.yaml` | bloqueado: reference-only canonico ausente |
| MLC SNP+Indel DITA windows | 2 | 0.96 | 0.96 | 0.96 | 0.96 | `snp_ancestry_predictor/configs/genotype_windows_all_variants_mle.yaml` | pronto, mas usa paths legados |
| CNN DITA AG only | 2 | 0.80 | 0.81 | 0.80 | 0.81 | `genotype_based_predictor/configs/genes_1000_all_3ontologies_signals_only.yaml` | pronto |

## Glossario Dos Sinais

- `MLC`: Maximum Likelihood Classification.
- `AG`: sinal AlphaGenome RNA-seq base a base.
- `AG0`: AlphaGenome zerado em posicoes sem variacao em nenhum individuo; no config atual, `alphagenome_signal_variant_mask: true`.
- `AGN`: AlphaGenome normalizado em relacao ao individuo de referencia; no config atual, `alphagenome_signal_transform: delta_reference`.
- `MI`: mascara de insercao.
- `MD`: mascara de delecao.
- `MV`: mascara de validade.
- `MS`: mascara de SNP.
- `DITA`: Dynamic INDEL Tensor Alignment usando o chain file do `bcftools consensus` para alinhar predicoes AlphaGenome ao eixo expandido global.

## Comandos De Reprodução

### Genotype CNN/PCA

Treinar:

```bash
python3 -m genomics_cli genotype train <config.yaml>
```

Avaliar checkpoint PyTorch:

```bash
python3 -m genomics_cli genotype evaluate <config.yaml> --checkpoint best_accuracy --split test
```

Os resultados ficam em:

```text
results/genotype_based_predictor/runs/<run_name>/
```

Arquivos principais:

```text
manifest.json
config.yaml
resolved_config.json
models/best_accuracy.pt
models/best_loss.pt
models/training_history.json
```

### MLC/SNP

Rodar:

```bash
cd snp_ancestry_predictor
python3 snp_ancestry_predictor.py --config configs/<config>.yaml
```

O pipeline MLC executa, conforme YAML:

- gerar arquivos 23andMe por individuo;
- computar frequencias alelicas no conjunto de referencia;
- predizer ancestralidade/target no conjunto de avaliacao.

## Mapeamento Detalhado

### MLC SNP Genome-Wide Panel

Resultado historico: `0.97` accuracy com 2 haplotipos.

Partir de:

```text
snp_ancestry_predictor/configs/default.yaml
```

Ajustes necessarios para repetir exatamente:

```yaml
statistics:
  snp_panel: null  # ou painel genome-wide usado originalmente
prediction:
  haplotype_mode: "H1+H2"
  method: "mle"
```

Status: precisa criar um config canonico versionado, porque o YAML atual mistura caminhos legados e painel de janelas 32k.

### MLC SNP Gene-Window Panel

Resultado historico: `0.71` accuracy com 2 haplotipos.

Partir de:

```text
snp_ancestry_predictor/configs/default.yaml
```

Ajustes:

```yaml
statistics:
  snp_panel: "/dados/GENOMICS_DATA/top3/refs/panel_23andme_V5_genes_1000_w32768.txt"
prediction:
  haplotype_mode: "H1+H2"
  method: "mle"
```

Status: precisa criar config canonico e migrar paths legados.

### PCA+RF

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies_pca400_rf.yaml
```

Notas:

- Usa `model.type: RF`.
- Usa PCA com `pca_components: 400`.
- Requer `scikit-learn`, `scipy` e `joblib`.

### PCA+XGBoost

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies_pca400_xgboost.yaml
```

Notas:

- Usa `model.type: XGBOOST`.
- Usa PCA com `pca_components: 400`.
- Requer `xgboost`, alem de `scikit-learn`, `scipy` e `joblib`.

### CNN Raw Center Crop

Baseline legado sem DITA:

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/neural_legacy/default.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/neural_legacy/default.yaml --checkpoint best_accuracy --split test
```

Para 2 haplotipos:

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/neural_legacy/genes_1000_all.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/neural_legacy/genes_1000_all.yaml --checkpoint best_accuracy --split test
```

Esses configs usam:

```yaml
tensor_layout: "raw_center_crop"
feature_mode: "signals_only"
```

### CNN DITA AG+MI+MD+MV

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy --split test
```

### CNN DITA MI+MD+MV

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies_masks_only.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/genes_1000_all_3ontologies_masks_only.yaml --checkpoint best_accuracy --split test
```

### CNN DITA MI+MD+MV+MS

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies_masks_snp.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/genes_1000_all_3ontologies_masks_snp.yaml --checkpoint best_accuracy --split test
```

### CNN DITA AG0+MI+MD+MV+MS

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies_variant_signal_mask.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/genes_1000_all_3ontologies_variant_signal_mask.yaml --checkpoint best_accuracy --split test
```

### CNN DITA AGN+MI+MD+MV+MS

Config:

```text
genotype_based_predictor/configs/genes_1000_all_3ontologies_delta_reference.yaml
```

Status atual: inativo.

Bloqueio:

```yaml
reference_predictions_dataset_dir: "/dados/GENOMICS_DATA/top3/reference_only_results_genes_1000_all_3ontologies"
```

Para reativar:

1. Materializar um dataset reference-only canonico.
2. Atualizar `reference_predictions_dataset_dir` para o novo path canonico.
3. Remover `metadata.active: false`.
4. Rodar `python3 -m genomics_cli audit-configs --fail-on-active-legacy`.

### MLC SNP+Indel DITA Windows

```bash
cd snp_ancestry_predictor
python3 snp_ancestry_predictor.py --config configs/genotype_windows_all_variants_mle.yaml
```

Config relevante:

```yaml
conversion:
  variant_types: ["snp", "indel"]
  genotype_encoding: "gt"
statistics:
  max_snps: null
  variant_types: ["snp", "indel"]
prediction:
  haplotype_mode: "H1+H2"
  method: "mle"
```

Status: reproduz a familia do resultado historico, mas ainda usa paths legados em `/dados/GENOMICS_DATA/top3`. Deve ser migrado para paths canonicos em uma etapa futura.

### CNN DITA AG Only

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/genes_1000_all_3ontologies_signals_only.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/genes_1000_all_3ontologies_signals_only.yaml --checkpoint best_accuracy --split test
```

## Template Para Novos Experimentos Importantes

Copie esta seção quando um experimento deixar de ser exploratorio e passar a ser resultado de referencia.

````markdown
### Nome curto do experimento

Data: YYYY-MM-DD
Branch/commit: `<git commit>`
Objetivo: ...
Dataset: ...
Config: `path/to/config.yaml`
Comando de treino:

```bash
...
```

Comando de avaliação:

```bash
...
```

Run dir:

```text
results/...
```

Resultado test:

| Precision | Recall | F1 | Accuracy |
| ---: | ---: | ---: | ---: |
|  |  |  |  |

Notas:

- ...
````

## Gaps De Reprodutibilidade

- Criar configs canonicos para os dois experimentos MLC SNP historicos `(1)` e `(2)` sem depender de edicao manual de `default.yaml`.
- Migrar configs `snp_ancestry_predictor` de `/dados/GENOMICS_DATA/top3` para dataset/VCFs canonicos ou registrar esses insumos no `genomics_pipeline.data_registry`.
- Materializar reference-only AlphaGenome predictions canonicas para reativar `AGN`.
- Registrar, para cada resultado importante novo, o commit exato e o diretório de run em `results/`.
