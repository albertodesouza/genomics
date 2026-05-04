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
