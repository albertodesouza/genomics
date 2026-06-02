# Research Pipelines

Este repositório agora tem uma CLI comum para executar etapas separadas sem decorar os entrypoints internos de cada pacote.

Resultados e comandos de experimentos importantes ficam em `docs/IMPORTANT_EXPERIMENTS.md`.

## Layout padrão

- `/dados/GENOMICS_DATA`: fonte padrão dos datasets brutos, materializados e processados, compartilhada entre máquinas.
- `results/`: resultados de experimentos, checkpoints, dashboards exportados e relatórios dentro deste repositório para inspeção pela IDE.
- `results/cache/`: caches intermediários reutilizáveis entre experimentos.

`results/` fica fora do Git.

O root dos datasets pode ser sobrescrito com `GENOMICS_DATA_ROOT`. O root de resultados pode ser sobrescrito com `GENOMICS_RESULTS_ROOT`, mas o default é `./results` no repositório.

## CLI comum

```bash
python3 -m genomics_cli --help
python3 -m genomics_cli audit-configs
python3 -m genomics_cli audit-data --dataset-id 1kg_high_coverage
```

## Estrutura de cada run

Os pipelines integrados à infraestrutura comum criam runs em `results/<pipeline>/runs/<run_name>/` com:

- `manifest.json`: metadados do run, ambiente, status e artefatos principais.
- `config.yaml`: cópia do YAML usado na execução.
- `resolved_config.json`: configuração validada/expandida usada pelo código.
- `models/`: checkpoints, histórico de treino e artefatos de modelo.
- `logs/`: reservado para logs de execução.
- `plots/`: reservado para figuras.
- `reports/`: reservado para relatórios.

Os pipelines ativos integrados a essa estrutura sao `genotype_based_predictor` e `variant_transformer_predictor`. O `neural_ancestry_predictor_deprecated` fica apenas como legado para reproducibilidade historica; novos experimentos devem ser migrados para `genotype_based_predictor/configs/neural_legacy/`.

## Componentes compartilhados

- `genomics_workspace.py`: resolve roots padrão de dados e resultados.
- `genomics_pipeline.data_registry`: resolve IDs lógicos de datasets para paths em `/dados/GENOMICS_DATA`.
- `genomics_pipeline.experiment`: cria estrutura de runs e `manifest.json`.
- `genomics_pipeline.config_io`: leitura YAML/JSON, escrita JSON e hashes estáveis.
- `genomics_pipeline.dataset_metadata`: leitura canônica de `dataset_metadata.json`, catálogo de janelas, pedigree e fontes VCF.
- `genomics_pipeline.reproducibility`: seeds Python/NumPy/PyTorch e `worker_init_fn`.
- `genomics_pipeline.splitting`: agrupamento family-aware e geração de splits compartilháveis.
- `genomics_pipeline.arrays`: utilitários NumPy compartilhados, como `block_reduce_2d`.
- `genomics_pipeline.sklearn_pca_cache`: cache compartilhado de `StandardScaler + PCA` para baselines sklearn.
- `genomics_pipeline.metrics`: metricas de classificacao, reports e JSONs de resultado.
- `genomics_pipeline.optim`: factory compartilhada de otimizadores.
- `genomics_pipeline.checkpointing`: save/load/resolve de checkpoints.
- `genomics_pipeline.targets`: mapeamento de targets diretos e derivados.
- `genomics_pipeline.run_utils`: selecao de device, split loader e campos comuns de manifest.
- `genomics_pipeline.data_loading`: kwargs padronizados de `DataLoader` e geradores deterministas.
- `genomics_pipeline.torch_collate`: padding tensorial para batches variaveis.
- `genomics_pipeline.training_utils`: `EpochTrainer` generico, schedulers e historico de treino compartilhados.
- `genomics_pipeline.wandb_utils`: inicializacao/finalizacao opcional de W&B.

`neural_ancestry_predictor_deprecated/sklearn_pca_cache.py` é apenas um wrapper temporário de compatibilidade. Código novo deve importar `genomics_pipeline.sklearn_pca_cache`.

## Fonte única dos dados

O dataset canônico para novas análises é `/dados/GENOMICS_DATA/v1/1kG_high_coverage`, exposto em `genomics_workspace.DEFAULT_DATASET_DIR` e `CANONICAL_1KG_HIGH_COVERAGE_DIR`.

Código novo pode referenciar o dataset por ID lógico em vez de path absoluto:

```yaml
dataset_input:
  dataset_id: "1kg_high_coverage"
```

ou, no caso de derivados:

```yaml
dataset:
  source_dataset_id: "1kg_high_coverage"
```

O registry reconhece:

- `1kg_high_coverage`: dataset canônico em `v1/1kG_high_coverage`.
- `legacy_top3_1kg_high_coverage`: dataset legado historico `top3/non_longevous_results_genes_1000_all`, mantido apenas como referencia de proveniencia.
- `variant_transformer_superpopulation`: derivado em `variant_transformer/superpopulation`.
- `variant_transformer_superpopulation_32k`: derivado em `variant_transformer/superpopulation_32k`.
- `variant_transformer_pigmentation_binary`: derivado em `variant_transformer/pigmentation_binary`.

O antigo diretório `/dados/GENOMICS_DATA/top3` foi movido para quarentena em `/dados/GENOMICS_DATA/_deprecated/top3_20260601`. Ele nao deve ser usado em novas analises.

`top3` está descontinuado para novas análises. Use o audit para localizar usos remanescentes:

```bash
python3 -m genomics_cli audit-configs --legacy-only
python3 -m genomics_cli audit-configs --fail-on-legacy
python3 -m genomics_cli audit-configs --fail-on-active-legacy
```

Use o audit fisico para validar datasets registrados e artefatos esperados:

```bash
python3 -m genomics_cli audit-data --dataset-id 1kg_high_coverage --fail-on-missing
python3 -m genomics_cli audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 3 --fail-on-missing
```

O segundo comando retorna código diferente de zero enquanto houver qualquer config com path `legacy-top3`. O terceiro comando falha apenas para configs ativos.

Status atual dos blockers de legado:

- `genotype_based_predictor` com `alignment_mapping: bcftools_chain` ja tem os artefatos necessarios materializados no dataset canonico por hardlink: `*.window.raw.fa`, `*.window.consensus_ready.vcf.gz`, indices `.tbi`, `*.window.vcf.gz` e indices `.tbi`.
- `genes_1000_all_3ontologies_delta_reference.yaml` continua inativo porque aponta `reference_predictions_dataset_dir` para um dataset legado de referencia ainda nao materializado no canonico.
- Alguns configs dentro de `neural_ancestry_predictor_deprecated/configs/` continuam apontando para paths historicos; isso e aceitavel apenas porque o pacote esta deprecated.

Critério mínimo para manter novas execucoes livres de legado: `python3 -m genomics_cli audit-configs --fail-on-active-legacy` deve passar. Esse audit passou apos a quarentena e a materializacao dos artefatos `bcftools_chain` no canonico.

Datasets derivados, como `/dados/GENOMICS_DATA/variant_transformer/*`, devem ser materializados a partir do dataset canônico e registrar a fonte em seus metadados.

### genotype_based_predictor

```bash
python3 -m genomics_cli genotype prepare-cache genotype_based_predictor/configs/repo_layout.example.yaml
python3 -m genomics_cli genotype train genotype_based_predictor/configs/repo_layout.example.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/repo_layout.example.yaml --checkpoint best_accuracy --split test
python3 -m genomics_cli genotype workbench
python3 -m genomics_cli genotype single-gene-screen genotype_based_predictor/configs/neural_legacy/pigmentation_binary_single_gene_screen.yaml --dry-run
```

Layouts de tensor principais:

- `haplotype_channels`: layout canonico alinhado para experimentos novos, normalmente com `alignment_mapping: bcftools_chain`.
- `raw_center_crop`: baseline legado sem alinhamento, usado pelos configs migrados de `neural_ancestry_predictor_deprecated` em `genotype_based_predictor/configs/neural_legacy/`.

### variant_transformer_predictor

```bash
python3 -m genomics_cli variant materialize --dataset-id 1kg_high_coverage --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation
python3 -m genomics_cli variant train variant_transformer_predictor/configs/repo_layout.example.yaml
python3 -m genomics_cli variant evaluate variant_transformer_predictor/configs/repo_layout.example.yaml --checkpoint best_accuracy --split test
python3 -m genomics_cli variant analyze-counts /dados/GENOMICS_DATA/variant_transformer/superpopulation
```

### neural_ancestry_predictor_deprecated

Este pipeline e legado/monolitico. A CLI comum apenas delega para o script preservado, e seu uso deve ficar restrito a reproducibilidade historica.

```bash
python3 -m genomics_cli neural train neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml
python3 -m genomics_cli neural test neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml
python3 -m genomics_cli neural pca-cache neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml
```

## Migração dos configs históricos

Os configs históricos podem continuar apontando para `/dados`. O que deve migrar para o repositório é o destino de resultados: `results_dir` para runs e `processed_cache_dir` quando o cache for artefato intermediário de experimento.

Durante a transição, os comandos da CLI aceitam overrides como `--dataset-dir`, `--processed-cache-dir`, `--results-dir`, `--processed-dir` e `--consensus-dataset-dir`, gerando uma cópia temporária do YAML para execução.

Overrides opcionais da CLI so sobrescrevem o YAML quando a flag e passada explicitamente. Isso permite que configs temporarios ou de smoke mantenham `results_dir` e caches proprios sem serem substituidos pelos defaults da CLI.

## Smokes validados

Smokes reais curtos foram executados para os dois pipelines ativos:

- `genotype_based_predictor`: treino de 1 epoca com `raw_center_crop`, `debug.max_samples_per_epoch: 16`, checkpoint `best_accuracy` e avaliacao em `test`.
- `variant_transformer_predictor`: treino de 1 epoca com modelo reduzido (`d_model=64`, `layers=1`, `heads=4`) e avaliacao em `test`.

Artefatos gerados:

- `results/genotype_based_predictor/runs_smoke/cnn2_pigmentation_rna_seq_H1_raw_center_crop_32768_log_s1k6x32f16_s2f32_s3f64_gpavg_fc256_L100-40_relu_0.5_adam/`
- `results/variant_transformer_predictor/runs_smoke/variant_transformer_superpopulation_d64_l1_h4_63deddaac4/`
