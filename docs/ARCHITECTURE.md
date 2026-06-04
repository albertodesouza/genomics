# Repository Architecture

Este repositorio contem mais de um pipeline. A organizacao alvo separa codigo por papel operacional, evitando que pipelines ativos, infraestrutura compartilhada, ferramentas auxiliares e codigo legado fiquem todos no root.

## Papeis Principais

```text
genomics
  core                  infraestrutura compartilhada
  workflows             pipelines de processamento e builders de dataset
  predictors            modelos e pipelines de predicao
  converters            conversores reutilizaveis
legacy                  codigo historico mantido para reproducibilidade
third_party             projetos externos vendorizados ou modificados
native                  ferramentas compiladas / nao Python
configs                 presets de execucao por dominio
docs                    documentacao transversal
scripts                 operacao, ambiente, diagnostico e manutencao
tests                   testes automatizados
```

## Layout Alvo

```text
src/
  genomics/
    __init__.py
    cli.py
    workspace.py

    core/
      config_io.py
      data_loading.py
      data_registry.py
      dataset_metadata.py
      experiment.py
      metrics.py
      optim.py
      reproducibility.py
      run_utils.py
      splitting.py
      targets.py
      checkpointing.py
      torch_collate.py
      torch_utils.py
      training_utils.py
      wandb_utils.py
      arrays.py
      sklearn_pca_cache.py

    workflows/
      genomes_analyzer/
      dataset_builders/
        non_longevous/
        longevity/
      alphagenome/

    predictors/
      genotype_based/
      variant_transformer/
      snp_ancestry/

    converters/
      vcf_to_23andme/

legacy/
  neural_ancestry_predictor_deprecated/
  neural_longevity_dataset/        # enquanto depender do layout top3 legado

third_party/
  FROGAncestryCalc/

native/
  genes_difference_count/
```

## Mapa Dos Modulos Atuais

| Modulo atual | Papel atual | Destino alvo |
|---|---|---|
| `src/genomics/core/` | Infraestrutura compartilhada | codigo ativo |
| `src/genomics/workspace.py` | Paths padrao e roots | codigo ativo |
| `src/genomics/cli.py` | CLI comum `genomics` | codigo ativo |
| `src/genomics/workflows/genomes_analyzer/` | Pipeline operacional FASTQ/BAM/CRAM/VCF | codigo ativo |
| `src/genomics/predictors/genotype_based/` | Predictor denso/alinhado | codigo ativo |
| `src/genomics/predictors/variant_transformer/` | Predictor de variantes esparsas | codigo ativo |
| `src/genomics/predictors/snp_ancestry/` | Predictor SNP/frequencias | codigo ativo |
| `src/genomics/converters/vcf_to_23andme/` | Conversor VCF para 23andMe | codigo ativo |
| `src/genomics/workflows/alphagenome/` | Integracao AlphaGenome | codigo ativo |
| `src/genomics/workflows/dataset_builders/non_longevous/` | Builder historico 1000G/AlphaGenome | codigo ativo |
| `legacy/neural_longevity_dataset/` | Builder de longevidade legado | legado |
| `legacy/neural_ancestry_predictor_deprecated/` | Predictor deprecated | legado |
| `third_party/FROGAncestryCalc/` | Projeto externo modificado | `third_party/FROGAncestryCalc/` |
| `native/genes_difference_count/` | Ferramenta C++/OpenMP | `native/genes_difference_count/` |

## Regras De Dependencia

- `genomics.core` nao deve importar `genomics.predictors`, `genomics.workflows` ou `genomics.converters`.
- `genomics.predictors` pode importar `genomics.core` e `genomics.workspace`.
- `genomics.workflows` pode importar `genomics.core`, `genomics.workspace` e `genomics.converters`.
- `genomics.converters` deve ser reutilizavel e evitar depender de predictors especificos.
- Codigo novo nao deve importar modulos em `legacy/`.
- Codigo novo nao deve depender diretamente de `third_party/`; crie adaptadores nos workflows quando necessario.

## Imports Alvo

```python
from genomics.core.metrics import classification_metrics
from genomics.workspace import DEFAULT_DATASET_DIR
from genomics.predictors.genotype_based.config import PipelineConfig
from genomics.predictors.variant_transformer.model import VariantTransformerClassifier
from genomics.converters.vcf_to_23andme.converter import convert_vcf_to_23andme
```

## CLI Alvo

```bash
genomics genotype train ...
genomics variant train ...
genomics genomes-analyzer run -c configs/genomes_analyzer/human_30x.yaml
genomics snp-ancestry run -c configs/predictors/snp_ancestry/default.yaml
genomics convert vcf-to-23andme -c configs/converters/vcf_to_23andme/default.yaml
```

Comandos antigos por wrappers Python foram removidos. Use o console script `genomics` ou `python3 -m genomics`.
