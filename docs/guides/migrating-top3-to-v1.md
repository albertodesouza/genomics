# Migrando o `top3` para o layout v1

Este guia converte os datasets do diretório antigo `/dados/GENOMICS_DATA/top3` para o layout canônico usado hoje pelo repositório (`/dados/GENOMICS_DATA/v1`), preservando também os datasets de genes aleatórios, e diz o que pode ser apagado para liberar espaço.

A conversão **move** os arquivos em vez de copiá-los: não usa hardlinks nem symlinks, e quando origem e destino estão no mesmo filesystem ela praticamente não ocupa espaço extra. A origem é consumida e no fim sobra só metadado, que se apaga.

## O que vai para onde

| Origem (`top3/`) | Destino | Conteúdo |
|---|---|---|
| `non_longevous_results_genes_1000_all` | `v1/1kG_high_coverage` | 11 genes de pigmentação × 3202 indivíduos (dataset canônico, `dataset_id: 1kg_high_coverage`) |
| `non_longevous_results_genes_1000_random` | `v1/1kG_high_coverage_random` | 11 genes aleatórios × 1300 indivíduos |
| `non_longevous_results_genes_1000_random_11_1` | `v1/1kG_high_coverage_random_11_1` | 11 genes aleatórios × 1300 indivíduos |
| `non_longevous_results_genes_1000_random_11_2` | `v1/1kG_high_coverage_random_11_2` | 11 genes aleatórios × 1300 indivíduos |
| `longevity_dataset/vcf_chromosomes/` | `v1/1kG_high_coverage/raw_variants/vcf_chromosomes/` | VCFs por cromossomo do 1000 Genomes high coverage |
| `refs/GRCh38_full_analysis_set_plus_decoy_hla.fa` (+ `.fai`) | `references/` | Referência GRCh38 |

Os três datasets aleatórios ficam em diretórios próprios, e não dentro de `1kG_high_coverage`, porque são outra coorte (1300 indivíduos) e foram gerados por outro build. Os configs de treino apontam para eles com `dataset_input.dataset_dir`.

O layout v1 de cada dataset fica assim:

```text
<dataset>/
  dataset_metadata.json
  layout_metadata.json
  gtf_cache.feather, selected_samples.csv, metadata_statistics.json, vcf_validation_report.json
  references/windows/<gene>/ref.window.fa
  references/windows/<gene>/window_metadata.json
  individuals/<amostra>/individual_metadata.json
  individuals/<amostra>/windows/<gene>/
    <amostra>.H1.window.fixed.fa, <amostra>.H2.window.fixed.fa
    <amostra>.H1.window.raw.fa, <amostra>.H2.window.raw.fa
    <amostra>.window.vcf.gz(.tbi), <amostra>.window.consensus_ready.vcf.gz(.tbi)
    predictions_H1/, predictions_H2/
```

A principal diferença para o formato antigo: o `ref.window.fa`, que antes se repetia dentro de cada indivíduo, passa a existir uma vez só por gene, em `references/windows/<gene>/`. Isso também libera espaço (~18 GB no dataset principal).

## 0. Preparação

Atualize o repositório e o ambiente. A opção `--move` do conversor é necessária:

```bash
cd <checkout do repositório genomics>
git pull
source scripts/env/start_genomics_universal.sh
python3 -m genomics.predictors.genotype_based.tools.materialize_dataset --help | grep -- --move
```

Defina as variáveis usadas no resto do guia:

```bash
OLD=/dados/GENOMICS_DATA/top3
V1=/dados/GENOMICS_DATA/v1
mkdir -p "$V1"
```

Confira se `top3` e `v1` estão no mesmo filesystem. As duas linhas devem mostrar o mesmo ponto de montagem:

```bash
df --output=target "$OLD" "$V1"
```

Se forem filesystems diferentes, a conversão ainda funciona: cada arquivo é copiado e apagado da origem antes de passar ao próximo, então o espaço extra é de um arquivo por vez. Só fica mais lenta.

## 1. O que pode apagar antes de converter

Nada nesta lista é usado pela conversão.

### Pode apagar

| Diretório | Tamanho | Por quê |
|---|---|---|
| `non_longevous_results_genes_1000` | 136G | Mesmos 11 genes de pigmentação, para os mesmos 1300 indivíduos dos datasets aleatórios. Tudo isso já está em `non_longevous_results_genes_1000_all`, com FASTAs e predições idênticas byte a byte (conferido por amostragem numa cópia do `top3`; confira na sua com o comando abaixo). |
| `non_longevous_results_genes` | 8G | 78 indivíduos com os mesmos 11 genes, também idênticos ao `genes_1000_all` |
| `non_longevous_results_dataset_cache` | 17M | Cache |
| `raw` | 33G | Arquivos `.sra` do `prefetch` do genomes-analyzer (trio NA12878), já extraídos para `fastq/` |
| `fastq_ds` | 22G | FASTQs subamostrados (intermediário do alinhamento) |
| `trimmed` | 27G | FASTQs após trimming (intermediário do alinhamento) |
| `refs/_bwa/*.tar.gz` | — | Tarball do índice BWA do genomes-analyzer. Só se o índice já estiver extraído (`refs/reference.fa.{amb,ann,bwt,pac,sa}` existem); com o índice presente, o pipeline não usa mais o tarball. |
| `caco_incompleto.tar.gz` | 1,3G | Só se o `caco.tar.gz` for a versão completa |

Antes de apagar, confira que os dois datasets são mesmo cópias do `genes_1000_all`. Faça isso **antes da conversão**, porque ela tira os arquivos do `genes_1000_all`. O comando compara 20 indivíduos sorteados em 4 genes e não deve imprimir nenhuma linha `DIFERENTE`:

```bash
for d in non_longevous_results_genes_1000 non_longevous_results_genes; do
  for i in $(ls "$OLD/$d/individuals" | shuf -n 20); do
    for g in TYR HERC2 MC1R OCA2; do
      for f in "$i.H1.window.fixed.fa" "$i.H2.window.fixed.fa" predictions_H1/rna_seq.npz; do
        cmp -s "$OLD/$d/individuals/$i/windows/$g/$f" "$OLD/non_longevous_results_genes_1000_all/individuals/$i/windows/$g/$f" \
          || echo "DIFERENTE: $d $i $g $f"
      done
    done
  done
done
```

Na raiz de um dataset (fora de `individuals/`) costumam ficar resultados de outras análises, como os de SNP ancestry. O conteúdo varia de máquina para máquina, então guarde **tudo o que estiver na raiz**, exceto `individuals/`, antes de apagar. É pouco espaço: metadados, JSONs e diretórios de resultados.

```bash
KEEP="$V1/1kG_high_coverage_runs/legacy_top3"
for d in non_longevous_results_genes_1000 non_longevous_results_genes; do
  mkdir -p "$KEEP/$d"
  find "$OLD/$d" -mindepth 1 -maxdepth 1 ! -name individuals -exec mv -t "$KEEP/$d/" {} +
  ls "$KEEP/$d"
done
```

Confira as listagens. Depois, apague:

```bash
rm -rf "$OLD/non_longevous_results_genes_1000" \
       "$OLD/non_longevous_results_genes" \
       "$OLD/non_longevous_results_dataset_cache" \
       "$OLD/raw" "$OLD/fastq_ds" "$OLD/trimmed"
rm -f  "$OLD"/refs/_bwa/*.tar.gz
```

### Também pode apagar, se não for realinhar

| Diretório | Tamanho | Observação |
|---|---|---|
| `fastq` | 70G | FASTQs do trio NA12878. Dá para baixar de novo do ENA (ERR3239334, ERR3989341, ERR3989342). Só são necessários para refazer o alinhamento, e o resultado dele está em `bam/`. |

### Decisão sua (não fazem parte do v1)

| Diretório | Tamanho | O que é |
|---|---|---|
| `non_longevous_results` | 78G | Primeiro experimento: 78 indivíduos, janelas de 1 Mb, ATAC + RNA-seq |
| `non_longevous_results_runs*` | ~130G | Treinos do preditor antigo (modelos, métricas, caches). **Mantenha os `_random*`** se quiser preservar esses experimentos. Nos outros, dá para guardar só métricas e figuras. |
| `longevity_dataset` (menos `vcf_chromosomes`) | ~60G | Dataset do pipeline de longevidade antigo (hoje em `legacy/`) |
| `bam`, `vcf`, `vep`, `fasta`, `genes`, `ancestry`, `comparisons`, `trio`, `paternity`, `qc`, `logs` | ~130G | Resultados do genomes-analyzer para o trio NA12878 |

## 2. Mover os VCFs e a referência

Faça isto antes da conversão: o conversor grava em cada janela o caminho do VCF do cromossomo e confere se ele existe.

Antes do `mv`, confira com `ls -l "$OLD/refs/reference.fa"` que o `reference.fa` do genomes-analyzer não é um symlink para `GRCh38_full_analysis_set_plus_decoy_hla.fa`. Se for, mover o arquivo quebra o genomes-analyzer: nesse caso, aponte o symlink para o caminho novo.

```bash
mkdir -p "$V1/1kG_high_coverage/raw_variants" /dados/GENOMICS_DATA/references
mv "$OLD/longevity_dataset/vcf_chromosomes" "$V1/1kG_high_coverage/raw_variants/"
mv "$OLD"/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   "$OLD"/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
   /dados/GENOMICS_DATA/references/
genomics references ensure-grch38
```

O `ensure-grch38` não baixa nada, porque o arquivo já está lá: ele só registra a referência (grava o `.metadata.json`).

## 3. Converter os quatro datasets

```bash
VCF="$V1/1kG_high_coverage/raw_variants/vcf_chromosomes/1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

python3 -m genomics.predictors.genotype_based.tools.materialize_dataset --move --vcf-pattern "$VCF" \
  "$OLD/non_longevous_results_genes_1000_all" "$V1/1kG_high_coverage"

for s in random random_11_1 random_11_2; do
  python3 -m genomics.predictors.genotype_based.tools.materialize_dataset --move --vcf-pattern "$VCF" \
    "$OLD/non_longevous_results_genes_1000_$s" "$V1/1kG_high_coverage_$s" || break
done
```

No fim, cada conversão imprime um resumo, por exemplo `(duplicates_removed=..., moved=...)`. O `duplicates_removed` conta as cópias repetidas de `ref.window.fa`.

- **Se aparecer `conflicts=N`:** N arquivos da origem eram diferentes do que já estava no destino e ficaram na origem. Não apague a origem antes de entender esses arquivos.
- **Se a conversão for interrompida** (queda, Ctrl-C): rode o mesmo comando de novo. Ela continua de onde parou.

## 4. Validar e apagar o que sobrou

Valide o dataset principal (confere todos os indivíduos):

```bash
genomics audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 0 --fail-on-missing
```

Em cada origem, os indivíduos devem ter ficado só com metadados e os marcadores `prediction_H*.ok.txt`. Os comandos abaixo não devem listar nada:

```bash
for d in non_longevous_results_genes_1000_all non_longevous_results_genes_1000_random \
         non_longevous_results_genes_1000_random_11_1 non_longevous_results_genes_1000_random_11_2; do
  find "$OLD/$d/individuals" -type f ! -name individual_metadata.json ! -name 'prediction_H?.ok.txt'
done
```

Na raiz de cada origem ainda podem estar resultados de outras análises, como os de SNP ancestry (`ancestry_results_*`, `snp_ancestry_statistics_*.json`) que os configs `configs/predictors/snp_ancestry/pigmentation_binary*.yaml` gravam dentro do dataset. A conversão leva para o v1 só os arquivos do próprio dataset (`gtf_cache.feather`, `selected_samples.csv`, `metadata_statistics.json`, `vcf_validation_report.json`). Por isso, **guarde o resto da raiz** antes de apagar:

```bash
KEEP="$V1/1kG_high_coverage_runs/legacy_top3"
for d in non_longevous_results_genes_1000_all non_longevous_results_genes_1000_random \
         non_longevous_results_genes_1000_random_11_1 non_longevous_results_genes_1000_random_11_2; do
  mkdir -p "$KEEP/$d"
  find "$OLD/$d" -mindepth 1 -maxdepth 1 ! -name individuals -exec mv -t "$KEEP/$d/" {} +
  ls "$KEEP/$d"
done
```

Se o `find` anterior não listou nada e as listagens acima estão certas, pode apagar:

```bash
rm -rf "$OLD/non_longevous_results_genes_1000_all" \
       "$OLD/non_longevous_results_genes_1000_random" \
       "$OLD/non_longevous_results_genes_1000_random_11_1" \
       "$OLD/non_longevous_results_genes_1000_random_11_2"
```

O resultado do `du -sh "$V1"/*` deve ficar próximo da soma dos quatro datasets antigos, menos as cópias repetidas de `ref.window.fa`.

## 5. Depois da migração

**Configs de treino.** O dataset principal é encontrado por `dataset_id: "1kg_high_coverage"`. Para os aleatórios, use o diretório. Os treinos antigos em `non_longevous_results_runs_genes_1000_random*` apontavam para `top3`; para reproduzi-los, troque o caminho:

```yaml
dataset_input:
  dataset_dir: "/dados/GENOMICS_DATA/v1/1kG_high_coverage_random_11_1"
```

**Configs do construtor de datasets.** Os YAMLs em `configs/workflows/non_longevous_dataset/` ainda usam os caminhos antigos dos VCFs e da referência. Isso só importa se você for construir janelas novas. Nesse caso, troque:

| Campo | Caminho antigo | Caminho novo |
|---|---|---|
| `vcf_pattern` | `/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes/...` | `/dados/GENOMICS_DATA/v1/1kG_high_coverage/raw_variants/vcf_chromosomes/...` |
| `fasta` | `/dados/GENOMICS_DATA/top3/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa` | `/dados/GENOMICS_DATA/references/GRCh38_full_analysis_set_plus_decoy_hla.fa` |

**Genes de controle do experimento de especificidade (opcional).** No v1 da máquina do Breno, o `1kG_high_coverage` tem 44 genes: os 11 de pigmentação, mais os 33 genes aleatórios reconstruídos para a coorte de pigmentação (1072 indivíduos). Para tê-los também, há duas opções:

- **Gerar de novo:** rodar `configs/workflows/non_longevous_dataset/specificity_control_genes*.yaml`, que precisa do AlphaGenome e leva horas de GPU.
- **Copiar do Breno:** fazer rsync das pastas `individuals/*/windows/<gene>` e `references/windows/<gene>` a partir da máquina dele.
