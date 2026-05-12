# Genotype Based Predictor

Modulo para explorar o dataset genomico local, construir views, visualizar tracks AlphaGenome, comparar sequencias alinhadas com INDELs e inspecionar experimentos de classificacao.

O uso recomendado hoje e pelo Workbench web local. Ele abre varias ferramentas em uma unica porta e evita que o usuario precise navegar por arquivos manualmente.

## 1. Ambiente

Antes de rodar comandos pesados ou viewers que consultam VCF, ative o ambiente `genomics`:

```bash
source scripts/start_genomics_universal.sh
```

Esse ambiente deve disponibilizar `bcftools`. Sem `bcftools`, alguns comandos ainda funcionam por fallback Python, mas consultas a `.vcf.gz` ficam muito mais lentas.

Dataset padrao usado pelas ferramentas:

```text
/dados/GENOMICS_DATA/v1/1kG_high_coverage
```

Dataset preservado com os FASTAs/VCFs usados para gerar as predicoes AlphaGenome:

```text
/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all
```

Esse diretorio contem `ref.window.fa`, `*.window.consensus_ready.vcf.gz`, `*.window.raw.fa`, `*.window.fixed.fa` e `predictions_H1/H2`. Ele e necessario para o modo padrao `bcftools_chain`, que reproduz o consenso do `bcftools` e gera chain files para mapear as predicoes AlphaGenome ao eixo expandido usado no treinamento.

Runs padrao:

```text
/dados/GENOMICS_DATA/v1/1kG_high_coverage_runs
```

TSVs alinhados padrao para o viewer:

```text
genotype_based_predictor/aligned_dna_genes_1000_all
```

## 2. Workbench Web

Inicie tudo com:

```bash
source scripts/start_genomics_universal.sh
python3 -m genotype_based_predictor.genomics_workbench
```

Abra no navegador:

```text
http://127.0.0.1:8780
```

Rotas principais:

```text
http://127.0.0.1:8780/apps/dataset/
http://127.0.0.1:8780/apps/view-builder/
http://127.0.0.1:8780/apps/experiments/
http://127.0.0.1:8780/apps/tracks/
http://127.0.0.1:8780/apps/alignment/
```

Com paths explicitos:

```bash
python3 -m genotype_based_predictor.genomics_workbench \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --runs-root /dados/GENOMICS_DATA/v1/1kG_high_coverage_runs \
  --aligned-tsv-root genotype_based_predictor/aligned_dna_genes_1000_all \
  --host 127.0.0.1 \
  --port 8780
```

Para acessar de outra maquina na rede local:

```bash
python3 -m genotype_based_predictor.genomics_workbench --host 0.0.0.0 --port 8780
```

Logs dos subprocessos ficam em:

```text
/tmp/genomics_workbench_logs
```

## 3. Ferramentas Do Workbench

### Dataset Browser

Arquivo:

```text
dataset_browser.py
```

Uso direto:

```bash
python3 -m genotype_based_predictor.dataset_browser \
  /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --host 127.0.0.1 \
  --port 8770
```

Serve para explorar individuos, genes, populacoes, superpopulacoes e metadados do dataset. Os filtros principais usam checkboxes com busca.

### View Builder

Arquivo:

```text
view_builder.py
```

Uso direto:

```bash
python3 -m genotype_based_predictor.view_builder \
  /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --host 127.0.0.1 \
  --port 8771
```

Serve para criar arquivos `.view.json` sem editar JSON manualmente. Permite selecionar genes, populacoes, superpopulacoes e individuos por interface visual.

### Experiment Dashboard

Arquivo:

```text
experiment_dashboard.py
```

Uso direto:

```bash
python3 -m genotype_based_predictor.experiment_dashboard \
  /dados/GENOMICS_DATA/v1/1kG_high_coverage_runs \
  --host 127.0.0.1 \
  --port 8772
```

Mostra experimentos, metricas, configs, resultados JSON e matrizes de confusao quando disponiveis.

### AlphaGenome Track Viewer

Arquivo:

```text
alphagenome_track_viewer.py
```

Uso direto:

```bash
python3 -m genotype_based_predictor.alphagenome_track_viewer \
  /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --host 127.0.0.1 \
  --port 8774
```

Permite comparar tracks AlphaGenome por gene, individuo, haplotipo, output e track. Quando `align with bcftools_chain` esta ativo, o viewer usa o mesmo mapeamento FASTA/predicao -> eixo expandido do treinamento.

### Aligned DNA Viewer

Arquivo:

```text
aligned_dna_viewer.py
```

Uso direto:

```bash
python3 -m genotype_based_predictor.aligned_dna_viewer \
  genotype_based_predictor/aligned_dna_genes_1000_all \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --alignment-mapping bcftools_chain \
  --consensus-dataset-dir /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all \
  --host 127.0.0.1 \
  --port 8765
```

O viewer le os TSVs gerados por `export_aligned_dna.py` para obter o eixo expandido, a referencia e metadados de janela. No modo padrao `bcftools_chain`, as sequencias dos individuos renderizadas na tabela sao reconstruidas a partir dos FASTAs reais usados pelo AlphaGenome e dos chain files gerados por `bcftools consensus -c`. Assim, a visualizacao usa o mesmo mapeamento que o treinamento.

Recursos atuais:

- selecao de gene;
- filtros por superpopulacao, populacao e individuo;
- selecao de `H1`, `H2` ou ambos;
- navegacao por janela;
- destaque de mutacoes em vermelho;
- destaque de `X`/gap em amarelo;
- coluna `index`, que e a posicao no eixo alinhado/expandido;
- coluna `ref genome pos`, que e a posicao no genoma de referencia;
- celula vazia em `ref genome pos` quando o `REF` tem `X`;
- tabela de variantes do VCF presentes na janela selecionada;
- destaque verde da faixa que sera enviada para a CNN durante treinamento;
- linha azul no centro do gene mapeado para o eixo expandido;
- botao `Ir para janela da CNN` para abrir diretamente a faixa de treino.

Limite padrao:

```text
200000 celulas por requisicao
```

Para alterar:

```bash
python3 -m genotype_based_predictor.aligned_dna_viewer \
  genotype_based_predictor/aligned_dna_genes_1000_all \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --alignment-mapping bcftools_chain \
  --consensus-dataset-dir /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all \
  --max-cells 300000
```

Evite selecionar muitos individuos com janelas grandes. Isso pode travar o navegador mesmo quando o backend esta protegido.

## 4. Mapeamento FASTA AlphaGenome Para Eixo Expandido

O treinamento usa `dataset_input.alignment_mapping`. O modo recomendado e:

```yaml
dataset_input:
  alignment_mapping: "bcftools_chain"
  consensus_dataset_dir: "/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all"
```

Nesse modo, para cada gene/individuo/haplotipo, o pipeline:

- roda `bcftools consensus -c` usando o `ref.window.fa` e o `*.window.consensus_ready.vcf.gz` preservados;
- valida que o consenso reconstruido bate com `*.window.raw.fa`;
- aplica a mesma regra de `adjust_to_target_size` e valida contra `*.window.fixed.fa`;
- parseia o chain para mapear indices do FASTA/predicao AlphaGenome para posicoes da referencia;
- usa o `DynamicIndelAligner` para definir o eixo expandido global e os slots de insercao;
- gera um `entry` com `copy_from_indices` absolutos no FASTA/predicao e `expanded_indices` no eixo global;
- copia os valores das predicoes AlphaGenome para o tensor final usando esse mapa.

O modo antigo `dynamic_indel` ainda existe para compatibilidade e diagnostico, mas nao deve ser usado como fonte principal quando as predicoes AlphaGenome vieram dos FASTAs gerados por `bcftools consensus`.

Cache do mapeador chain:

```text
/dados/GENOMICS_DATA/v1/1kG_high_coverage/alignment_cache/bcftools_chain_mapper_v3
```

## 5. Gerar TSVs De DNA Alinhado

O exportador usa o `DynamicIndelAligner` para projetar REF, H1 e H2 no mesmo eixo expandido. As bases de H1/H2 no TSV sao reconstruidas diretamente do VCF cromossomico e da referencia. Esse TSV e util para inspecionar o eixo e as variantes, mas no modo recomendado `bcftools_chain` as sequencias dos individuos usadas para treino/visualizacao sao derivadas dos FASTAs reais do AlphaGenome e dos chain files.

Isso e importante porque os FASTAs fixos foram gerados para alimentar o AlphaGenome com comprimento constante. Eles podem ter sido truncados/padronizados apos `bcftools consensus` e nao devem ser usados como fonte de verdade para validar coordenadas base-a-base contra o VCF.

Arquivo:

```text
export_aligned_dna.py
```

Exemplo pequeno:

```bash
source scripts/start_genomics_universal.sh
python3 -m genotype_based_predictor.export_aligned_dna \
  genotype_based_predictor/configs/genes_1000_all.yaml \
  /tmp/aligned_dna_MC1R_check.tsv \
  --gene MC1R \
  --sample-limit 5 \
  --batch-size 2
```

Por padrao, o export alinha apenas a regiao central usada pelo modelo, com `--center-window-size 32768`. Isso evita consultar todos os VCFs da janela de aproximadamente 512 kbp quando o objetivo e validar/treinar a CNN de 32k.

Para voltar ao comportamento antigo e alinhar a janela completa:

```bash
python3 -m genotype_based_predictor.export_aligned_dna \
  genotype_based_predictor/configs/genes_1000_all.yaml \
  /tmp/aligned_dna_MC1R_full.tsv \
  --gene MC1R \
  --sample-limit 5 \
  --full-window
```

Formato gerado:

```text
# gene=MC1R
# expanded_length=524906
sample_id    H1_aligned    H2_aligned
REF          ...           ...
HG00096      ...           ...
```

O caractere `X` significa ausencia de base naquela coluna alinhada para aquele haplotipo. Em uma delecao `CA -> C`, por exemplo, a base ancora `C` fica preservada e o `X` fica na linha da base removida `A`.

### Export Em Lotes

Para todos os individuos, sempre use `--batch-size`. Isso limita memoria.

Recomendacao para maquina com 16 GB RAM:

```text
--batch-size 2 ou --batch-size 4
```

Exemplo seguro para um gene:

```bash
source scripts/start_genomics_universal.sh
mkdir -p genotype_based_predictor/aligned_dna_genes_1000_all

nice -n 10 python3 -m genotype_based_predictor.export_aligned_dna \
  genotype_based_predictor/configs/genes_1000_all.yaml \
  genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.tsv \
  --gene MC1R \
  --all-samples \
  --batch-size 4
```

Exemplo para os 11 genes da view `genes_1000_all`:

```bash
source scripts/start_genomics_universal.sh
mkdir -p genotype_based_predictor/aligned_dna_genes_1000_all

for gene in MC1R TYRP1 TYR SLC45A2 DDB1 EDAR MFSD12 OCA2 HERC2 SLC24A5 TCHH; do
  nice -n 10 python3 -m genotype_based_predictor.export_aligned_dna \
    genotype_based_predictor/configs/genes_1000_all.yaml \
    "genotype_based_predictor/aligned_dna_genes_1000_all/${gene}.tsv" \
    --gene "$gene" \
    --all-samples \
    --batch-size 4
done
```

Se o computador ficar lento, interrompa com `Ctrl+C` e reduza para:

```text
--batch-size 2
```

## 6. Converter TSV Para Texto Em Colunas

Arquivo:

```text
aligned_dna_columns.py
```

Uso:

```bash
python3 -m genotype_based_predictor.aligned_dna_columns \
  genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.tsv \
  genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.color.columns.txt \
  --color
```

Abrir preservando cores ANSI:

```bash
less -R genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.color.columns.txt
```

Para muitos individuos, prefira o `Aligned DNA Viewer`. O arquivo em colunas pode ficar gigantesco.

## 7. DynamicIndelAligner

Arquivo:

```text
dynamic_indel_alignment.py
```

Responsabilidades:

- consultar variantes de janela no VCF;
- construir eixo de referencia expandido;
- preservar coordenadas homologas entre individuos;
- mapear delecoes para `X` na base removida, preservando a base ancora do VCF;
- mapear insercoes para colunas extras apos a base ancora;
- expor `get_alignment_axis(...)` e `get_haplotype_entry(...)`.

Detalhe importante: o aligner valida o `REF` do VCF contra `ref.window.fa`. Quando a posicao bruta do VCF nao bate exatamente com a sequencia local, ele procura um pequeno deslocamento local e usa a posicao que corresponde ao `REF`. Isso evita desalinhamentos como o caso `CA -> C`, em que o `X` deve cair na linha do `A`, nao na linha do `C` ancora.

### Janela Central Da CNN

Para treino CNN, a janela de `window_center_size` nao e cortada pelo centro geometrico do eixo expandido. A politica correta e:

1. localizar o centro da janela do gene na referencia original;
2. converter esse indice de referencia para o eixo expandido;
3. cortar uma janela simetrica no eixo expandido ao redor desse ponto.

Assim, para `window_center_size=32768`, o aligner consulta somente essa regiao central da referencia por padrao. O eixo expandido resultante pode ter mais que 32.768 colunas quando ha insercoes dentro da regiao, mas continua centrado no gene. Isso evita que um numero assimetrico de INDELs antes/depois do centro desloque a regiao biologica analisada.

Essa politica e registrada na cache processada como:

```text
center_window_policy = reference_center_to_expanded_axis
```

Caches processadas antigas sem essa politica sao invalidadas.

### Cache Compartilhada De Alinhamento

O `DynamicIndelAligner` e a fonte unica de verdade para alinhamento usado por:

- treino/CNN;
- `export_aligned_dna.py`;
- `AlphaGenome Track Viewer` quando alinhamento dinamico esta ativo;
- geracao de TSVs que o `Aligned DNA Viewer` consome.

A cache persistente padrao fica dentro do dataset:

```text
/dados/GENOMICS_DATA/v1/1kG_high_coverage/alignment_cache/dynamic_indel_ref_window_v2/
```

A cache e separada por conjunto de amostras, porque o eixo expandido depende das insercoes observadas no conjunto analisado. O namespace tem o formato:

```text
samples_<N>_<hash>/<GENE>/
```

Exemplo da cache do eixo expandido global:

```text
alignment_cache/dynamic_indel_ref_window_v4/samples_5_d8dd50e1ef446666/MC1R/
```

Arquivos principais:

```text
axis.json
samples/HG00096.json
samples/HG00097.json
```

O `axis.json` guarda o eixo expandido do gene. Cada arquivo em `samples/` guarda as entradas H1/H2 de um individuo naquele mesmo eixo.

No modo recomendado `bcftools_chain`, existe tambem a cache do mapeamento FASTA AlphaGenome -> eixo expandido:

```text
alignment_cache/bcftools_chain_mapper_v3/<GENE>/<AXIS_KEY>/<SAMPLE>/<HAP>.entry.json
```

Essa entrada e gerada a partir de `bcftools consensus -c`, validada contra os FASTAs existentes e usada tanto no treinamento quanto na visualizacao.

A cache inclui metadados de proveniencia com:

- versao do algoritmo de alinhamento;
- dataset usado;
- conjunto de amostras;
- `window_metadata.json`;
- `ref.window.fa`;
- VCF usado.

Se esses metadados mudarem, a entrada antiga deixa de ser usada.

### Relação Com A Cache Processada Da CNN

A cache processada da CNN continua em:

```text
<processed_cache_dir>/datasets/<dataset_name>/
```

mas agora ela salva em `metadata.json` uma `alignment_cache_signature`. Caches processadas antigas, sem essa assinatura, sao consideradas invalidas e recriadas.

Isso garante que o treino nao reutilize tensores gerados por um algoritmo de alinhamento antigo.

## 8. Views E Configs

Views ficam em:

```text
genotype_based_predictor/views/
```

Configs ficam em:

```text
genotype_based_predictor/configs/
```

Regra pratica:

- `views/`: define o que entra no modelo ou nos viewers, como genes, individuos, populacoes e outputs.
- `configs/`: define como o experimento usa esses dados, como modelo, treino, split e metrica.

Campos importantes para alinhamento das predicoes AlphaGenome:

```yaml
dataset_input:
  alignment_mapping: "bcftools_chain"
  consensus_dataset_dir: "/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all"
```

Use `dynamic_indel` apenas para diagnostico/compatibilidade com caches antigas.

Config atual de referencia para os 11 genes principais:

```text
genotype_based_predictor/configs/genes_1000_all.yaml
```

View associada:

```text
genotype_based_predictor/views/genes_1000_all.view.json
```

Genes:

```text
MC1R TYRP1 TYR SLC45A2 DDB1 EDAR MFSD12 OCA2 HERC2 SLC24A5 TCHH
```

## 9. Treino E Experimentos

O treinamento principal continua em:

```text
train.py
```

Antes de treinar, valide a entrada do pipeline `bcftools_chain`:

```bash
source scripts/start_genomics_universal.sh
python3 -m genotype_based_predictor.validate_training_input \
  genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml \
  --sample-limit 5 \
  --max-tensor-items 3
```

O validador reconstrói consensos com `bcftools consensus -c`, valida os FASTAs preservados, checa limites dos `.npz` AlphaGenome e monta alguns tensores finais pelo `ProcessedGenomicDataset`. Para a config de 3 ontologias, o shape esperado deve ser fixo, por exemplo:

```text
(2, 66, 32768)
```

Para acelerar a primeira execução, precompute a cache `bcftools_chain` em paralelo antes do treino:

```bash
source scripts/start_genomics_universal.sh
python3 -m genotype_based_predictor.precompute_bcftools_chain_cache \
  genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml \
  --workers 6 \
  --chunk-size 25
```

Esse comando e retomavel: entradas ja geradas sao reaproveitadas. Se uma reconstrução cacheada nao bater com `*.window.raw.fa` ou `*.window.fixed.fa`, o mapper regenera o consenso/chain uma vez antes de falhar.

Se quiser processar em blocos retomáveis:

```bash
python3 -m genotype_based_predictor.precompute_bcftools_chain_cache \
  genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml \
  --workers 6 \
  --start 0 \
  --limit 500
```

Exemplo geral:

```bash
source scripts/start_genomics_universal.sh
python3 -m genotype_based_predictor.train genotype_based_predictor/configs/genes_1000_all.yaml
```

Para a config de 3 ontologias:

```bash
source scripts/start_genomics_universal.sh
python3 -m genotype_based_predictor.train genotype_based_predictor/configs/genes_1000_all_3ontologies.yaml
```

Use o `Experiment Dashboard` para inspecionar resultados ja gerados.

## 10. Cuidados De Memoria

Regras praticas:

- Nao gere todos os genes em paralelo.
- Para export all-samples, use `--batch-size 2` ou `--batch-size 4` em maquinas com 16 GB RAM.
- Use `nice -n 10` para deixar o sistema mais responsivo durante exports longos.
- No viewer, evite selecionar centenas de individuos com janelas grandes.
- Se estiver usando o Workbench, prefira navegar por janelas pequenas e filtros.
- Regenerar TSVs invalida os `.idx.json` automaticamente se tamanho ou mtime mudarem.

## 11. Documentos Relacionados

- `docs/genomics_workbench_plan.md`: arquitetura e plano do Workbench.
- `docs/tratamento_indels.md`: explicacao tecnica sobre INDELs, eixo expandido e visualizacao.
- `views/README.md`: papel das views.
