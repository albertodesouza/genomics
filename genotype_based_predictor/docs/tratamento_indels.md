# Tratamento de Inserções e Deleções (INDELs) nas Saídas do AlphaGenome

## A. Visão Geral

No presente trabalho, cada indivíduo é representado a partir da saída do modelo AlphaGenome, obtida ao se fornecer como entrada uma sequência genômica haplotípica (H1 ou H2) construída em uma janela de 512 kbp centrada em um gene de interesse. A saída do modelo consiste em múltiplas trilhas funcionais ao longo da sequência, das quais extraímos uma janela central de 32.768 posições para cada gene.

Considerando 11 genes e 6 trilhas por gene, obtém-se uma matriz de dimensão $66 \times 32768$, utilizada como entrada para uma rede neural convolucional (CNN).

Entretanto, a presença de inserções e deleções (INDELs) entre indivíduos introduz desalinhamentos entre as posições correspondentes dessas matrizes, uma vez que as coordenadas das saídas do modelo passam a depender da sequência específica de cada indivíduo.

Para mitigar esse problema, propomos um método de alinhamento das saídas baseado em uma referência expandida, complementado por máscaras explícitas de inserção e deleção. O objetivo é permitir que posições homólogas entre indivíduos sejam comparadas de forma consistente, mesmo quando as respectivas sequências haplotípicas contêm INDELs.

---

## B. Exemplo Ilustrativo do Problema

Considere o seguinte exemplo simplificado, no qual uma pequena sequência de referência é comparada com três indivíduos:

```text
referência:    A C T A G T
individuo1:    A C G T           (2 deleções)
individuo2:    A C T T A G T     (1 inserção)
individuo3:    A C T A G G T     (1 inserção)
```

Nesse cenário, as sequências possuem comprimentos distintos e não podem ser comparadas posição a posição de forma direta. Em particular, inserções deslocam as posições subsequentes, enquanto deleções removem posições existentes na referência.

---

## C. Construção do Eixo de Referência Expandido

Definimos um eixo de referência expandido $R^*$, que acomoda todas as inserções observadas no conjunto analisado. No exemplo anterior, a referência expandida é dada por:

```text
referência*:   A C T - A G - T
```

Esse eixo possui 8 posições, em vez das 6 posições da referência original, incorporando:

- uma posição extra para representar a inserção observada no indivíduo 2;
- uma posição extra para representar a inserção observada no indivíduo 3.

As posições marcadas com `-` na referência expandida não correspondem a bases existentes na referência original; elas são espaços adicionais introduzidos apenas para permitir que inserções observadas em indivíduos sejam representadas em um eixo comum.

---

## D. Alinhamento das Sequências no Eixo Expandido

Cada indivíduo é então projetado nesse eixo comum:

```text
referência*:   A C T - A G - T
individuo1:    A C - - - G - T
individuo2:    A C T T A G - T
individuo3:    A C T - A G G T
```

Esse alinhamento garante que posições homólogas fiquem verticalmente alinhadas. No exemplo:

- o indivíduo 1 possui deleções nas posições correspondentes a `T` e `A` da referência original;
- o indivíduo 2 possui uma inserção na quarta coluna do eixo expandido;
- o indivíduo 3 possui uma inserção na sétima coluna do eixo expandido.

---

## E. Construção das Máscaras de INDELs

Para cada indivíduo, definimos máscaras binárias que indicam explicitamente a presença de inserções, deleções e posições válidas. Essas máscaras são usadas como canais adicionais de entrada para a rede neural.

### E.1 Indivíduo 1: deleções

```text
referência*:   A C T - A G - T
individuo1:    A C - - - G - T
máscara_i1:    0 0 0 0 0 0 0 0
máscara_d1:    0 0 1 0 1 0 0 0
máscara_val1:  1 1 0 0 0 1 0 1
```

Neste caso:

- `máscara_i1` é zero em todas as posições, pois o indivíduo 1 não possui inserções;
- `máscara_d1` é igual a 1 nas posições deletadas em relação à referência;
- `máscara_val1` indica onde existe valor real associado ao indivíduo.

### E.2 Indivíduo 2: inserção

```text
referência*:   A C T - A G - T
individuo2:    A C T T A G - T
máscara_i2:    0 0 0 1 0 0 0 0
máscara_d2:    0 0 0 0 0 0 0 0
máscara_val2:  1 1 1 1 1 1 0 1
```

Neste caso:

- `máscara_i2` é igual a 1 na posição da inserção;
- `máscara_d2` é zero em todas as posições;
- `máscara_val2` indica que há valor real em todas as posições ocupadas pelo indivíduo, exceto na coluna de expansão não usada por ele.

### E.3 Indivíduo 3: inserção em outra posição

```text
referência*:   A C T - A G - T
individuo3:    A C T - A G G T
máscara_i3:    0 0 0 0 0 0 1 0
máscara_d3:    0 0 0 0 0 0 0 0
máscara_val3:  1 1 1 0 1 1 1 1
```

Neste caso:

- `máscara_i3` é igual a 1 na posição da inserção;
- `máscara_d3` é zero em todas as posições;
- `máscara_val3` indica quais posições possuem valores reais para o indivíduo 3.

---

## F. Generalização para Saídas do AlphaGenome

No caso real, não trabalhamos diretamente com letras, mas sim com tensores contínuos produzidos pelo AlphaGenome.

Para cada indivíduo $i$, a saída original do AlphaGenome pode ser representada como:

$$
X_i \in \mathbb{R}^{C \times L_i}
$$

onde $C = 66$ corresponde aos canais resultantes da concatenação dos 11 genes e 6 tracks por gene, e $L_i$ é o comprimento efetivo da representação para o indivíduo, após considerar sua sequência haplotípica.

Utilizando os arquivos VCF e a referência GRCh38, construímos uma função de mapeamento entre coordenadas da referência e coordenadas do indivíduo:

$$
f_i : \mathrm{pos}_{ref} \rightarrow \mathrm{pos}_{ind}
$$

Esse mapeamento permite identificar, para cada posição do eixo expandido, se a posição corresponde a:

1. uma correspondência direta entre referência e indivíduo;
2. uma deleção no indivíduo;
3. uma inserção presente no indivíduo.

---

## G. Projeção das Saídas no Eixo Expandido

Cada saída $X_i$ é projetada no eixo expandido, produzindo um tensor alinhado:

$$
\tilde{X}_i \in \mathbb{R}^{C \times L^*}
$$

onde $L^*$ é o comprimento do eixo de referência expandido.

Para cada posição $j$ do eixo expandido:

- se a posição existe no indivíduo, copia-se para $\tilde{X}_i[:,j]$ o valor correspondente da saída do AlphaGenome;
- se a posição corresponde a uma deleção, atribui-se um valor neutro, por exemplo zero;
- se a posição corresponde a uma inserção do indivíduo, copia-se o valor da saída do AlphaGenome correspondente à base inserida;
- se a posição corresponde a uma inserção observada em outro indivíduo, mas ausente no indivíduo atual, atribui-se valor neutro.

Assim, o símbolo `-` do exemplo didático corresponde, no caso real, a uma posição do tensor preenchida por valor neutro e identificada por máscaras apropriadas.

---

## H. Representação Final para a CNN

A entrada final para a CNN é construída pela concatenação do tensor alinhado com as máscaras de inserção e deleção:

$$
Z_i = \mathrm{concat}(\tilde{X}_i,\ M_i^{(ins)},\ M_i^{(del)})
$$

resultando em:

$$
Z_i \in \mathbb{R}^{(C+2) \times L^*}
$$

onde:

- $\tilde{X}_i$ contém os valores alinhados derivados do AlphaGenome;
- $M_i^{(ins)}$ indica as posições de inserção;
- $M_i^{(del)}$ indica as posições de deleção.

Opcionalmente, pode-se também concatenar uma máscara de validade:

$$
Z_i = \mathrm{concat}(\tilde{X}_i,\ M_i^{(ins)},\ M_i^{(del)},\ M_i^{(val)})
$$

nesse caso resultando em:

$$
Z_i \in \mathbb{R}^{(C+3) \times L^*}
$$

A máscara de validade é útil para informar explicitamente à rede quais posições contêm valores reais e quais posições foram preenchidas por valores neutros.

---

## I. Exemplo Numérico Abstrato com Saídas do AlphaGenome

Para ilustrar a transição do exemplo com letras para o caso real, considere que cada letra represente um valor escalar produzido pelo AlphaGenome em uma determinada track. Nesse caso, poderíamos ter:

```text
referência*:       A   C   T   -   A   G   -   T
individuo2:        A   C   T   T   A   G   -   T
saída_track2:     x1  x2  x3  x4  x5  x6   0  x7
máscara_i2:        0   0   0   1   0   0   0   0
máscara_d2:        0   0   0   0   0   0   0   0
máscara_val2:      1   1   1   1   1   1   0   1
```

Para múltiplas tracks, o mesmo processo é realizado canal a canal. Assim, no lugar de um único valor $x_j$, a posição $j$ contém um vetor de dimensão $C$.

---

## J. Tratamento de Coordenadas e Bordas

Todas as operações devem ser realizadas em coordenadas genômicas absolutas. A janela de 32.768 posições deve ser definida em relação ao genoma de referência, e apenas INDELs que intersectam essa janela devem ser considerados.

Esse cuidado é especialmente importante em regiões próximas às bordas da janela. Por exemplo, uma inserção imediatamente antes do início da janela de interesse pode deslocar os índices locais da sequência individual, mas não deve deslocar artificialmente o eixo de análise definido sobre a referência.

Assim, a janela de interesse deve ser tratada como um intervalo genômico fixo, por exemplo:

```text
chr:start-end
```

e não apenas como um intervalo de índices locais dentro da sequência individual fornecida ao AlphaGenome.

---

## K. Inserções no Início da Janela

Um caso particular ocorre quando um ou mais indivíduos possuem inserções no início da janela central de 32.768 posições. Nessa situação, a referência expandida deve acomodar essas inserções antes da primeira posição de referência incluída na janela, desde que elas intersectem o intervalo considerado.

Por exemplo, se a janela de interesse começa na posição $p$ da referência e um indivíduo possui uma inserção imediatamente após $p$, o eixo expandido deve conter uma coluna adicional logo após essa posição. Indivíduos sem essa inserção recebem valor neutro nessa coluna e máscara de validade igual a zero.

Esse procedimento evita que a inserção no início da janela cause deslocamento global em todas as posições subsequentes.

---

## L. Extensão para Haplótipos H1 e H2

O procedimento descrito é aplicado separadamente aos dois haplótipos de cada indivíduo, H1 e H2. Existem pelo menos três formas naturais de incorporar os dois haplótipos na CNN:

1. concatenar H1 e H2 ao longo do eixo de canais;
2. processar H1 e H2 como entradas paralelas;
3. agregar as representações dos dois haplótipos antes da etapa classificatória.

Na primeira opção, se cada haplótipo produz uma representação de dimensão $(C+2) \times L^*$, a concatenação de H1 e H2 ao longo dos canais produz uma entrada de dimensão:

$$
2(C+2) \times L^*
$$

ou, caso a máscara de validade também seja usada:

$$
2(C+3) \times L^*
$$

---

## M. Generalização para Dados de Teste

O eixo expandido $R^*$ deve ser definido de modo compatível com o protocolo experimental. Uma escolha metodologicamente segura é construir $R^*$ usando apenas os conjuntos de treinamento e validação, evitando vazamento de informação do conjunto de teste.

Durante o teste, novos indivíduos são projetados nesse mesmo eixo. Caso surjam inserções não observadas durante o treinamento, há três possibilidades:

1. truncar essas inserções;
2. agregá-las à posição de referência vizinha;
3. usar uma estratégia dinâmica de extensão do eixo, desde que isso seja feito sem alterar a arquitetura já treinada ou sem invalidar a comparação experimental.

A primeira alternativa é a mais simples e preserva uma dimensão fixa de entrada. A terceira alternativa é mais flexível, mas exige cuidado adicional para manter a compatibilidade com a CNN.

---

## N. Discussão

A abordagem proposta transforma um problema de sequências de comprimento variável em uma representação tensorial alinhada, complementada por máscaras explícitas de INDELs. Dessa forma, diferenças estruturais entre indivíduos são preservadas como informação de entrada, em vez de serem eliminadas ou confundidas com deslocamentos artificiais.

Essa estratégia permite integrar modelos baseados em sequência, como o AlphaGenome, com arquiteturas discriminativas de aprendizado profundo, como CNNs ou transformers, preservando alinhamento posicional e coerência semântica entre indivíduos.

Além disso, o uso de máscaras torna explícito para o modelo quais posições correspondem a valores reais, inserções, deleções ou preenchimentos neutros. Isso reduz o risco de que a rede interprete valores artificiais como sinais biológicos genuínos.

---

## O. Resumo Operacional do Pipeline

O procedimento completo pode ser resumido da seguinte forma:

1. Para cada indivíduo, haplótipo e gene, construir a sequência haplotípica a partir do genoma de referência e dos arquivos VCF.
2. Executar o AlphaGenome normalmente sobre a sequência haplotípica.
3. Extrair a janela central de interesse.
4. Usar o VCF para construir o mapeamento entre coordenadas da referência e coordenadas do indivíduo.
5. Construir o eixo de referência expandido.
6. Projetar a saída do AlphaGenome no eixo expandido.
7. Preencher posições ausentes com valor neutro.
8. Construir as máscaras de inserção, deleção e, opcionalmente, validade.
9. Concatenar valores e máscaras para formar a entrada final da CNN.
10. Repetir o processo separadamente para H1 e H2.

### O.1 Normalizacao Local Do REF Do VCF

Na implementacao atual, o `DynamicIndelAligner` usa o arquivo `ref.window.fa` como fonte local de verdade para posicionar INDELs.

Para cada variante de comprimento diferente, o aligner compara o campo `REF` do VCF com a sequencia da janela. Se a posicao bruta do VCF nao bate exatamente com o `REF` local, ele procura um pequeno deslocamento local e usa a coordenada que realmente corresponde ao `REF`.

Esse passo evita deslocamentos artificiais em casos como:

```text
VCF: CA -> C
```

Nessa delecao, o `C` e a base ancora e deve permanecer alinhado. O `X` deve aparecer na linha da base removida `A`:

```text
REF:        C A
individuo: C X
```

O resultado correto no viewer fica equivalente a:

```text
540  C  C
541  A  X
```

e nao:

```text
540  C  X
541  A  C
```

Esse comportamento esta implementado em `dynamic_indel_alignment.py`.

### O.2 Cache Compartilhada De Alinhamento

O alinhamento dinamico e implementado uma unica vez em `DynamicIndelAligner` e reutilizado por treino, exportacao de TSVs e viewers.

A cache persistente padrao fica em:

```text
<dataset_dir>/alignment_cache/dynamic_indel_ref_window_v2/
```

Ela e separada por conjunto de amostras, pois o eixo expandido depende das insercoes observadas naquele conjunto. Portanto, uma exportacao com 5 individuos nao reaproveita o mesmo eixo de uma execucao com todos os individuos.

Dentro de cada conjunto, a estrutura e:

```text
samples_<N>_<hash>/<GENE>/axis.json
samples_<N>_<hash>/<GENE>/samples/<sample_id>.json
```

O `axis.json` guarda o eixo expandido. Os arquivos por amostra guardam `copy_from_indices`, `expanded_indices`, `insertion_indices` e `deletion_indices` para H1/H2.

A cache processada da CNN salva uma assinatura dessa cache de alinhamento em `metadata.json`. Caches processadas antigas sem essa assinatura sao invalidadas automaticamente.

### O.3 Janela 32k Centrada No Gene

Ao treinar a CNN, a janela central de 32.768 posicoes deve permanecer centrada no gene de interesse, definido no eixo da referencia, e nao no centro geometrico do eixo expandido.

A implementacao atual usa a seguinte regra:

1. define o centro da janela do gene no eixo de referencia original;
2. converte esse indice para o eixo expandido usando `expanded_index_map`;
3. corta uma janela simetrica ao redor desse indice expandido.

Para `window_center_size=32768`, isso produz 16.384 posicoes a esquerda e 16.384 posicoes a direita do centro de referencia mapeado para o eixo expandido, exceto se a janela bater em uma borda extrema.

Essa regra evita o erro de usar simplesmente `expanded_length // 2`, que poderia deslocar o centro biologico quando ha assimetria no numero de insercoes antes e depois do gene.

Por padrao, o `DynamicIndelAligner` consulta e alinha apenas essa regiao central de 32.768 bases da referencia. O comportamento antigo de consultar toda a janela de aproximadamente 512 kbp continua disponivel no export com `--full-window`.

Mesmo quando a consulta e limitada a 32.768 bases de referencia, o eixo expandido pode ficar ligeiramente maior que 32.768 colunas por causa de insercoes dentro da regiao.

A cache processada registra essa politica como:

```text
center_window_policy = reference_center_to_expanded_axis
```

---

## P. Visualização em Texto do Alinhamento de DNA

Para inspecionar se o `DynamicIndelAligner` está construindo o eixo expandido corretamente, existem dois utilitários no módulo:

```text
genotype_based_predictor/export_aligned_dna.py
genotype_based_predictor/aligned_dna_columns.py
```

O primeiro gera um arquivo TSV com a sequência alinhada da referência e dos indivíduos. O segundo converte esse TSV para uma visualização em colunas, em que cada posição alinhada ocupa uma linha e cada indivíduo/haplótipo ocupa uma coluna.

### P.1 Ambiente

Antes de executar os comandos, inicialize o ambiente `genomics`, que contém `bcftools`:

```bash
source scripts/start_genomics_universal.sh
```

Também é possível prefixar cada comando com `source scripts/start_genomics_universal.sh && ...`.

### P.2 Arquivo TSV alinhado

Exemplo para o gene `MC1R` e os 5 primeiros indivíduos da view `one_gene_10_individuals`:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.export_aligned_dna genotype_based_predictor/configs/one_gene_10_individuals.yaml genotype_based_predictor/aligned_dna_MC1R_ref_plus_5.tsv --sample-limit 5
```

O formato gerado é:

```text
# gene=MC1R
# expanded_length=524906
sample_id    H1_aligned    H2_aligned
REF          ...           ...
HG00096      ...           ...
HG00097      ...           ...
```

O caractere `X` indica uma coluna do eixo expandido que não possui base naquele haplótipo. Isso ocorre, por exemplo, quando a coluna representa uma inserção presente em outro indivíduo ou uma posição deletada no indivíduo atual.

### P.3 Visualização em colunas

Para converter o TSV para uma visualização comparável com `more` ou `less`:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.aligned_dna_columns genotype_based_predictor/aligned_dna_MC1R_ref_plus_5.tsv genotype_based_predictor/aligned_dna_MC1R_ref_plus_5.columns.txt
```

O formato fica:

```text
# gene=MC1R
# expanded_length=524906
pos    REF    HG00096_H1    HG00096_H2    HG00097_H1    HG00097_H2
1      G      G             G             G             G
2      C      C             C             C             C
3      C      C             C             C             C
```

### P.4 Visualização com cores

Para destacar diferenças em relação à referência e posições `X`, use `--color`:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.aligned_dna_columns genotype_based_predictor/aligned_dna_MC1R_ref_plus_5.tsv genotype_based_predictor/aligned_dna_MC1R_ref_plus_5.color.columns.txt --color
```

Convenção de cores:

- `X` aparece em amarelo.
- bases diferentes da referência aparecem em vermelho.
- bases iguais à referência ficam sem destaque.

Para abrir preservando as cores ANSI:

```bash
less -R genotype_based_predictor/aligned_dna_MC1R_ref_plus_5.color.columns.txt
```

### P.5 Todos os indivíduos de uma view

Para usar todos os indivíduos definidos pela view ou, se a view não define amostras explicitamente, todos os indivíduos de `dataset_metadata.json`, use `--all-samples`:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.export_aligned_dna genotype_based_predictor/configs/one_gene_10_individuals.yaml genotype_based_predictor/aligned_dna_MC1R_all.tsv --gene MC1R --all-samples --batch-size 4
```

O argumento `--batch-size` controla quantos individuos sao processados por vez. Em maquinas com 16 GB RAM, recomenda-se `--batch-size 2` ou `--batch-size 4` para evitar uso excessivo de memoria.

### P.6 Todos os genes de `genes_1000_all.yaml`

A configuração `genotype_based_predictor/configs/genes_1000_all.yaml` usa a view `genes_1000_all`, que inclui os genes:

```text
MC1R TYRP1 TYR SLC45A2 DDB1 EDAR MFSD12 OCA2 HERC2 SLC24A5 TCHH
```

Para gerar um TSV por gene, para todos os individuos do dataset, use execucao em lotes. Evite rodar varios genes em paralelo em maquinas com pouca memoria:

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

Se o computador ficar lento, interrompa e reduza para `--batch-size 2`.

Para gerar uma visualizacao colorida em texto de um gene ja exportado:

```bash
python3 -m genotype_based_predictor.aligned_dna_columns \
  genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.tsv \
  genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.color.columns.txt \
  --color
```

Para abrir um gene específico:

```bash
less -R genotype_based_predictor/aligned_dna_genes_1000_all/MC1R.color.columns.txt
```

### P.7 Observações de desempenho

O `DynamicIndelAligner` usa `bcftools query` quando `bcftools` está disponível no ambiente. Isso é o caminho recomendado para views grandes.

Se `bcftools` não estiver disponível, o código possui um fallback que lê o `.vcf.gz` diretamente em Python. Esse fallback é útil para inspeções pequenas, mas é mais lento e não é recomendado para todos os indivíduos ou muitos genes.

As visualizações em colunas podem ficar muito grandes. Para cerca de 1000 indivíduos, cada indivíduo gera duas colunas (`H1` e `H2`), além da coluna `REF` e da coluna `pos`. Portanto, cada arquivo por gene pode ter milhares de colunas e centenas de milhares de linhas. Nesses casos, prefira `less -R` em vez de `more`, ou gere subconjuntos menores com `--sample-limit`.

### P.8 Interface web para TSVs alinhados

Para explorar os TSVs sem converter tudo para uma tabela de texto gigante, use a interface web:

```text
genotype_based_predictor/aligned_dna_viewer.py
```

Ela lê arquivos `.tsv` gerados por `export_aligned_dna.py`, cria um índice pequeno (`.idx.json`) com os offsets das linhas e carrega apenas a janela de posições solicitada. Isso evita carregar o arquivo inteiro em memória.

Para iniciar a interface usando um diretório com um `.tsv` por gene:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.aligned_dna_viewer genotype_based_predictor/aligned_dna_genes_1000_all --host 127.0.0.1 --port 8765
```

Depois abra no navegador:

```text
http://127.0.0.1:8765
```

Funcionalidades:

- seleção do gene de interesse;
- seleção de indivíduos por checkboxes;
- filtros por superpopulação, população e busca textual;
- seleção de `H1`, `H2` ou ambos;
- escolha da posição inicial e tamanho da janela;
- navegação por janela anterior/próxima;
- opção de mostrar apenas posições com diferença em relação à referência;
- coluna `index`, correspondente ao eixo alinhado/expandido;
- coluna `ref genome pos`, correspondente à coordenada no genoma de referência;
- célula vazia em `ref genome pos` quando a linha do `REF` contém `X`;
- tabela de variantes do VCF contidas na janela selecionada;
- destaque visual de diferenças em vermelho;
- destaque de `X` em amarelo.

Por padrão, a API limita a renderização a 200.000 células por requisição. Isso protege máquinas com pouca memória e evita travar o navegador. Para alterar o limite:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.aligned_dna_viewer genotype_based_predictor/aligned_dna_genes_1000_all --host 127.0.0.1 --port 8765 --max-cells 300000
```

Evite selecionar todos os indivíduos com janelas muito grandes. Para uma análise ampla, use janelas menores ou ative a opção de mostrar apenas variantes.
