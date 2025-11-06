# Modificações no FROGAncestryCalc

## Resumo

O código fonte foi modificado para aceitar arquivos de entrada delimitados por **pipe (|)** em vez de vírgula (,).

## Arquivos Modificados

### 1. `src/dv/ValidateFileHeader.java`
- Alterado `split("\\,")` → `split("\\|")`
- Alterado `indexOf(',')` → `indexOf('|')`
- Alterado `countOccurrences(record, ',')` → `countOccurrences(record, '|')`
- Mensagens de erro atualizadas: "comma delimited" → "pipe delimited"

### 2. `src/read/ReadTxtFiles.java`
- Alterado `split("\\,")` → `split("\\|")` na leitura de genótipos

## Formato dos Arquivos de Entrada

Os arquivos de entrada devem ter o seguinte formato:

```
Individual|rs10497191|rs1079597|rs11652805|...(outros SNPs)
HG02561_GWD|NN|CC|CC|...(genótipos)
HG02562_GWD|TT|CT|CC|...(genótipos)
```

**Importante:**
- Delimitador: pipe `|` (não vírgula)
- Finais de linha: Unix (LF) 
- Encoding: UTF-8

## Como Executar

### A partir dos fontes compilados:

```bash
./run_from_source.sh
```

Ou manualmente:

```bash
LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 java -cp bin main.ComputeBatchAnalysis
```

### Nota sobre Locale

É essencial usar `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8` para evitar problemas com formatação de números em notação científica (ex: `5.652E-62` em vez de `5,652E-62`).

## Recompilação

Se precisar recompilar após modificações:

```bash
cd /home/lume2/genomics/frog/FROGAncestryCalc
rm -rf bin
mkdir bin
javac -d bin -sourcepath src $(find src -name "*.java")
cp -r src/read/data bin/read/
```

## Arquivos de Saída

A aplicação gera 3 arquivos em `output/`:
- `*_likelihood.txt` - Probabilidades de ancestralidade para cada população
- `*_orderOfMag.txt` - Ordem de magnitude das probabilidades
- `*_rankOrder.txt` - Ranking das populações por probabilidade

## Backup

Os arquivos originais foram salvos com extensão `.bak`:
- `src/dv/ValidateFileHeader.java.bak`

