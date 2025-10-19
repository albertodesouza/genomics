# 🧬 Como Baixar Sequências Genômicas Reais

## 📖 Visão Geral

Para usar o Neural Module com dados reais, você precisa baixar sequências de genomas de referência. Este guia mostra como fazer isso.

---

## 🌐 Opção 1: Ensembl REST API (Recomendado)

### Vantagens
- ✅ Rápido e simples
- ✅ Sem necessidade de instalar ferramentas
- ✅ Sempre atualizado

### Baixar Sequência por Coordenadas

```bash
# Formato: chr:start..end
# IMPORTANTE: Intervalos são inclusivos (start e end incluídos)
# Para 2048 bp: end = start + 2047

# Exemplo: Gene HBB (chr11:5227002-5229049 = 2048 bp)
curl 'https://rest.ensembl.org/sequence/region/human/11:5227002..5229049?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > HBB_gene.fasta

# Verificar tamanho
grep -v "^>" HBB_gene.fasta | tr -d '\n' | wc -c
# Deve mostrar: 2048
```

### Baixar por Nome de Gene

```bash
# Exemplo: Gene BRCA1
curl 'https://rest.ensembl.org/sequence/id/ENSG00000012048?type=genomic;coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > BRCA1_gene.fasta
```

### Assemblies Disponíveis
- **Humano**: `GRCh38` (hg38, atual) ou `GRCh37` (hg19, legado)
- **Camundongo**: `GRCm39` ou `GRCm38`
- **Outros**: Veja https://rest.ensembl.org/info/species

---

## 🔬 Opção 2: UCSC Genome Browser

### Via Interface Web

1. Acesse: https://genome.ucsc.edu/
2. Navegue para: **Tools → Table Browser**
3. Configure:
   - Assembly: `Human GRCh38/hg38`
   - Group: `Genes and Gene Predictions`
   - Track: `GENCODE v45`
   - Output format: `sequence`
4. Clique em **get output**

### Via Linha de Comando (com twoBitToFa)

```bash
# Instalar ferramenta (Ubuntu 20.04, ver http://hgdownload.cse.ucsc.edu/admin/exe/ para outros sistemas)
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v479/twoBitToFa
chmod +x twoBitToFa

# Baixar sequência
./twoBitToFa http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit:chr11:5227002-5229049 HBB_gene.fasta
```

---

## 🧪 Opção 3: NCBI (GenBank)

### Via Interface Web

1. Acesse: https://www.ncbi.nlm.nih.gov/gene/
2. Busque pelo gene (ex: "HBB human")
3. Na página do gene, vá em **Genomic regions, transcripts, and products**
4. Clique no link do GenBank (ex: NG_000007.3)
5. Clique em **Send to → File → FASTA**

### Via E-utilities

```bash
# Buscar gene HBB
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=HBB[Gene]+AND+human[Organism]&retmode=json"

# Baixar sequência por ID
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_000007.3&rettype=fasta&retmode=text" > HBB_gene.fasta
```

---

## 🔧 Opção 4: samtools + Genoma de Referência

Se você já tem o genoma de referência local:

```bash
# Baixar genoma de referência (uma vez)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

# Extrair região específica
samtools faidx hg38.fa chr11:5227002-5229049 > HBB_gene.fasta
```

---

## 📊 Ajustando para Tamanhos Suportados

O AlphaGenome suporta apenas tamanhos específicos:
- **2,048 bp** (2 KB)
- **16,384 bp** (16 KB)
- **131,072 bp** (128 KB)
- **524,288 bp** (512 KB)
- **1,048,576 bp** (1 MB)

### Ajustar Tamanho da Região

```bash
# Para 2048 bp: end = start + 2047
# Exemplo: chr1:1000000-1002047

# Para 16384 bp: end = start + 16383
# Exemplo: chr1:1000000-1016383
```

### Script de Validação

```bash
#!/bin/bash
# validate_sequence_size.sh

FASTA_FILE=$1

# Contar nucleotídeos (sem header)
SIZE=$(grep -v "^>" "$FASTA_FILE" | tr -d '\n' | wc -c)

# Tamanhos suportados
SUPPORTED=(2048 16384 131072 524288 1048576)

echo "Tamanho: $SIZE bp"

if [[ " ${SUPPORTED[@]} " =~ " ${SIZE} " ]]; then
    echo "✓ Tamanho suportado!"
else
    echo "✗ Tamanho não suportado!"
    echo "Tamanhos válidos: ${SUPPORTED[@]}"
fi
```

---

## 🎯 Exemplos Práticos

### Genes Humanos Comuns

```bash
# Gene HBB (Beta-globina) - Anemia Falciforme
curl 'https://rest.ensembl.org/sequence/region/human/11:5227002..5229049?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > HBB_2048bp.fasta

# Gene CFTR (Fibrose Cística) - 2048 bp da região promotora
curl 'https://rest.ensembl.org/sequence/region/human/7:117559590..117561637?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > CFTR_promoter_2048bp.fasta

# Gene TP53 (p53, supressor de tumor) - 2048 bp
curl 'https://rest.ensembl.org/sequence/region/human/17:7676154..7678201?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > TP53_2048bp.fasta

# Gene BRCA1 (Câncer de mama) - 16KB
curl 'https://rest.ensembl.org/sequence/region/human/17:43044295..43060678?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > BRCA1_16kb.fasta
```

### Regiões Regulatórias

```bash
# Enhancer no cromossomo 11 - 2048 bp
curl 'https://rest.ensembl.org/sequence/region/human/11:5200000..5202047?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > enhancer_chr11_2048bp.fasta
```

---

## ✅ Checklist de Uso

Antes de usar uma sequência baixada:

- [ ] Verificar tamanho: `grep -v "^>" seq.fasta | tr -d '\n' | wc -c`
- [ ] Confirmar que é um dos tamanhos suportados
- [ ] Verificar assembly correto (GRCh38 recomendado)
- [ ] Header do FASTA é informativo
- [ ] Sequência contém apenas A, T, C, G, N

---

## 🚀 Uso com Neural Module

```bash
# Baixar sequência
curl 'https://rest.ensembl.org/sequence/region/human/11:5227002..5229049?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > my_sequence.fasta

# Analisar
python neural_module.py \
  -i my_sequence.fasta \
  -k YOUR_API_KEY \
  -o results/

# Analisar variante
python neural_module.py \
  -i my_sequence.fasta \
  -k YOUR_API_KEY \
  -o results_variant/ \
  --variant 1024 A T

# Analisar com metadados de ontologia
python neural_module.py \
  -i my_sequence.fasta \
  -k YOUR_API_KEY \
  -o results_with_metadata/ \
  --save-metadata
```

---

## 🔗 Recursos Adicionais

- **Ensembl REST API Docs**: https://rest.ensembl.org/
- **UCSC Genome Browser**: https://genome.ucsc.edu/
- **NCBI Gene**: https://www.ncbi.nlm.nih.gov/gene/
- **AlphaGenome Docs**: https://www.alphagenomedocs.com/

---

**Criado**: Outubro 2025  
**Parte do**: Neural Module Documentation

