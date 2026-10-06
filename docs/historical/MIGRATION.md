# Repository Reorganization Migration Plan

Este documento descreve uma migracao em fases para sair do layout atual, com varios pacotes no root, para um pacote unico `genomics` com subdominios claros.

## Principios

- Fazer movimentos pequenos e verificaveis.
- Usar wrappers temporarios apenas durante movimentos intermediarios; eles devem ser removidos antes do layout final.
- Evitar alterar comportamento junto com mudancas de layout.
- Atualizar testes e documentacao na mesma fase do movimento correspondente.
- Remover compatibilidade antiga somente quando nao houver consumidores externos ou scripts dependentes.

## Fase 1: Documentacao E Classificacao

Status: concluida nesta etapa.

Objetivos:

- Documentar a arquitetura alvo em `docs/historical/ARCHITECTURE.md`.
- Atualizar o `README.md` raiz para explicar o repo como plataforma multi-pipeline.
- Manter o codigo no lugar por enquanto.

Validacao:

```bash
python3 -m pytest tests
```

## Fase 2: Criar Pacote `genomics`

Status: concluida. As implementacoes reais ficam em `src/genomics/cli.py` e `src/genomics/workspace.py`; wrappers antigos foram removidos.

Objetivos:

- Adicionar `pyproject.toml`.
- Criar `src/genomics/__init__.py`.
- Mover a implementacao de `genomics_cli.py` para `src/genomics/cli.py`.
- Mover a implementacao de `genomics_workspace.py` para `src/genomics/workspace.py`.
- Criar `src/genomics/core/` como namespace canonico para a infraestrutura antes mantida em `genomics_pipeline/`.
- Criar `genomics/` no root como shim para permitir `python3 -m genomics` e imports `genomics.*` direto do checkout sem instalacao.
- Migrar imports dos predictors ativos e testes de `genomics_pipeline.*` para `genomics.core.*`.
- Remover wrappers no root depois da consolidacao. Concluido.

Imports alvo:

```python
from genomics.workspace import DEFAULT_DATASET_DIR
```

Validacao:

```bash
python3 -m genomics.cli --help
python3 -m genomics --help
python3 -m pytest tests/test_genomics_cli.py tests/test_data_registry.py
```

## Fase 3: Mover Infra Compartilhada Para `genomics.core`

Status: concluida. Os arquivos reais foram movidos para `src/genomics/core/`; wrappers antigos foram removidos.

Objetivos:

- Mover `genomics_pipeline/*` para `src/genomics/core/`. Concluido.
- Atualizar imports internos dos pipelines ativos. Concluido para `genomics.predictors.genotype_based` e `genomics.predictors.variant_transformer`.
- Remover wrappers temporarios em `genomics_pipeline/`. Concluido.

Imports alvo:

```python
from genomics.core.metrics import classification_metrics
from genomics.core.training_utils import EpochTrainer
```

Validacao:

```bash
python3 -m pytest tests/test_training_utils.py tests/test_experiment_run.py tests/test_model_forward_smoke.py
```

## Fase 4: Mover Predictors Ativos

Status: concluida para os predictors ativos Python. `variant_transformer_predictor/` foi movido para `src/genomics/predictors/variant_transformer/`. `genotype_based_predictor/` foi movido para `src/genomics/predictors/genotype_based/`. Wrappers antigos foram removidos; configs canonicos ficam em `configs/predictors/`.

Objetivos:

- Mover `genotype_based_predictor/` para `src/genomics/predictors/genotype_based/`. Concluido; wrappers de compatibilidade removidos.
- Mover `variant_transformer_predictor/` para `src/genomics/predictors/variant_transformer/`. Concluido.
- Atualizar `genomics.cli` para chamar os novos modulos. Concluido para `variant_transformer` e `genotype_based`.
- Remover wrappers temporarios para os entrypoints antigos. Concluido.

Comandos alvo:

```bash
genomics genotype train configs/predictors/genotype_based/repo_layout.example.yaml
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
```

Validacao:

```bash
python3 -m pytest tests/test_genotype_package_layout.py tests/test_model_forward_smoke.py
genomics genotype --help
genomics variant --help
```

## Fase 5: Mover Workflows, Converters E Ferramentas

Status: concluida para codigo Python. Workflows, converters e predictors auxiliares foram movidos para `src/genomics/`; wrappers Python antigos foram removidos. `genes_difference_count/` foi movido para `native/genes_difference_count/`. `FROGAncestryCalc/` foi movido para `third_party/FROGAncestryCalc/`.

Objetivos:

- Mover `genomes_analyzer_pipeline/` para `src/genomics/workflows/genomes_analyzer/`. Concluido.
- Mover `neural_module/` para `src/genomics/workflows/alphagenome/`. Concluido para codigo Python.
- Mover `build_non_longevous_dataset/` para `src/genomics/workflows/dataset_builders/non_longevous/`. Concluido para codigo Python.
- Mover `snp_ancestry_predictor/` para `src/genomics/predictors/snp_ancestry/`. Concluido para codigo Python.
- Mover `vcf_to_23andme/` para `src/genomics/converters/vcf_to_23andme/`. Concluido.
- Mover `genes_difference_count/` para `native/genes_difference_count/`. Concluido.
- Mover `FROGAncestryCalc/` para `third_party/FROGAncestryCalc/`. Concluido.

Validacao:

```bash
genomics genomes-analyzer run --help
genomics dataset-builders non-longevous build --help
genomics alphagenome analyze -- --help
genomics convert vcf-to-23andme --help
python3 -m pytest tests
```

## Fase 6: Isolar Legado

Status: concluido para codigo Python. Implementacoes reais foram movidas para `legacy/`; wrappers Python antigos foram removidos.

Objetivos:

- Mover `neural_ancestry_predictor_deprecated/` para `legacy/neural_ancestry_predictor_deprecated/`. Concluido para codigo Python.
- Mover `neural_longevity_dataset/` para `legacy/neural_longevity_dataset/` enquanto depender de `/dados/GENOMICS_DATA/top3`. Concluido para codigo Python.
- Atualizar docs para indicar que codigo novo nao deve importar esses modulos.

Validacao:

```bash
genomics audit-configs --fail-on-active-legacy
```

## Fase 7: Reorganizar Configs, Scripts E Docs

Status: em andamento. A arvore canonica `configs/` foi criada por dominio e recebeu copias dos YAMLs historicos. Os arquivos originais permanecem nos modulos antigos como assets de compatibilidade ate a migracao de scripts externos e documentos historicos. Scripts foram migrados para `scripts/env/`, `scripts/ops/`, `scripts/maintenance/`, `scripts/diagnostics/`, `scripts/experiments/` e `scripts/dev/`; os caminhos historicos em `scripts/` permanecem como wrappers de compatibilidade.

Objetivos:

- Reorganizar `configs/` por dominio.
- Reorganizar `scripts/` em `env/`, `ops/`, `maintenance/`, `diagnostics/`, `experiments/` e `dev/`.
- Espelhar `tests/` de acordo com `src/genomics/`.
- Atualizar todos os READMEs com os comandos novos.

Layout alvo de configs:

```text
configs/
  genomes_analyzer/
  predictors/
    genotype_based/
    variant_transformer/
    snp_ancestry/
  workflows/
    non_longevous_dataset/
    longevity_dataset/
    alphagenome/
  converters/
    vcf_to_23andme/
```

## Criterio De Conclusao

A migracao pode ser considerada concluida quando:

- Codigo novo usa apenas imports `genomics.*`.
- Wrappers antigos foram removidos.
- `python3 -m pytest tests` passa.
- `docs/historical/ARCHITECTURE.md` reflete o layout real.
- O README raiz e um indice da plataforma, nao a documentacao de um unico pipeline.
