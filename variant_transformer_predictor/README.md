# Variant Transformer Predictor

Transformer esparso para classificacao baseada apenas em variantes geneticas de VCF faseado. O modulo nao usa AlphaGenome, nao cria representacao densa do genoma e nao alinha genomas.

## Fluxo

1. Materializar tokens esparsos a partir de VCFs em `/dados`.
2. Treinar o Transformer com RoPE baseado em `position_relative` real.
3. Avaliar checkpoints em train/val/test.

## Materializacao

Usando regioes/metadados de um dataset ja materializado do repo:

```bash
source scripts/start_genomics_universal.sh
python3 -m variant_transformer_predictor.materialize_dataset \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --vcf-root-dir /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes \
  --target superpopulation \
  --classes AFR AMR EAS EUR SAS \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation \
  --sample-batch-size 128
```

Para materializar apenas a janela central de 32 kb de cada gene/janela:

```bash
python3 -m variant_transformer_predictor.materialize_dataset \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --vcf-root-dir /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes \
  --target superpopulation \
  --classes AFR AMR EAS EUR SAS \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation_32k \
  --sample-batch-size 128 \
  --central-window-size 32768
```

Para analisar a quantidade de tokens por gene antes/depois do recorte central:

```bash
python3 -m variant_transformer_predictor.analyze_variant_counts \
  /dados/GENOMICS_DATA/variant_transformer/superpopulation \
  --central-window-size 32768
```

Tambem e possivel usar BED e metadados TSV/JSON explicitamente:

```bash
python3 -m variant_transformer_predictor.materialize_dataset \
  --regions-bed /dados/GENOMICS_DATA/regions/pigmentation_genes.bed \
  --samples-metadata /dados/GENOMICS_DATA/metadata.tsv \
  --vcf-pattern '/dados/GENOMICS_DATA/vcf/1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz' \
  --target superpopulation \
  --classes AFR AMR EAS EUR SAS \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation
```

Formato BED esperado:

```text
chrom  start0  end0  gene_id
```

O BED e interpretado como 0-based half-open. O VCF e consultado como 1-based inclusive.

## Treino

```bash
python3 -m variant_transformer_predictor.train variant_transformer_predictor/configs/superpopulation.yaml
```

## Avaliacao

```bash
python3 -m variant_transformer_predictor.evaluate_checkpoint \
  variant_transformer_predictor/configs/superpopulation.yaml \
  --checkpoint best_accuracy \
  --split test
```

## Representacao

Cada variante observada em um haplotipo gera um token com:

```text
variant_type, position_relative, haplotype, gene, ref_allele, alt_allele, length_norm
```

O modelo adiciona `[CLS]` internamente e aplica RoPE somente em `Q` e `K` da atencao, usando `position_relative`.

## Observacoes

- VCFs nao faseados sao ignorados por padrao (`unphased_policy=skip`).
- `max_sequence_length` controla o limite de tokens por individuo. Para genome-wide, o custo de atencao `O(N^2)` pode ficar alto.
- A materializacao consulta VCF por regiao e por lotes de amostras via `bcftools query`; sem `bcftools`, usa fallback gzip mais lento.
