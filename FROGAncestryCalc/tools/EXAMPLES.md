# Exemplos de Uso das Ferramentas

Este arquivo contém exemplos práticos de como usar as ferramentas de extração de SNPs.

## Exemplo 1: Converter arquivo existente do 1000 Genomes

Se você já tem um arquivo VCF do 1000 Genomes:

```bash
# Supondo que você baixou manualmente um VCF
python3 tools/vcf_to_frog.py \
    /path/to/1000genomes.vcf.gz \
    tools/aisnps_55_list.txt \
    input/my_1000g_samples.txt

# Configurar para análise
echo "inputFilename=my_1000g_samples.txt" >> FROGAncestryCalc.properties
echo "panelInfo=55AI" >> FROGAncestryCalc.properties

# Executar análise
./run.sh
```

## Exemplo 2: Baixar dados do 1000 Genomes automaticamente

```bash
# Baixar todos os dados (usa VCFs existentes ou baixa se necessário)
./tools/extract_snps_from_1000genomes.sh

# Ou baixar apenas algumas amostras específicas
cat > test_samples.txt << 'SAMPLES'
HG02561
HG02562
HG03055
SAMPLES

./tools/extract_snps_from_1000genomes.sh \
    -s test_samples.txt \
    -o input/test_samples.txt

# Analisar
./run.sh
```

## Exemplo 3: Processar dados de sequenciamento clínico

### De VCF (mais comum)

```bash
# Assumindo que você tem sample.vcf.gz do sequenciamento
./tools/extract_snps_from_wgs.sh \
    -i /path/to/patient_001.vcf.gz \
    -t vcf \
    -o input/patient_001.txt

# Atualizar configuração
sed -i 's/inputFilename=.*/inputFilename=patient_001.txt/' FROGAncestryCalc.properties

# Executar análise
./run.sh
```

### De BAM (se você só tem arquivo alinhado)

```bash
# Baixar referência (se não tiver)
# wget http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

./tools/extract_snps_from_wgs.sh \
    -i /path/to/sample.bam \
    -t bam \
    -r /path/to/GRCh38.fa \
    -o input/sample.txt \
    -n MySample

# Executar análise
./run.sh
```

### De FASTQ (sequenciamento bruto)

```bash
./tools/extract_snps_from_wgs.sh \
    -i sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -t fastq \
    -r /path/to/GRCh38.fa \
    -o input/sample.txt \
    -n MySample
```

## Exemplo 4: Processar múltiplas amostras em lote

```bash
# Se você tem vários VCFs na mesma pasta
for vcf in vcf_files/*.vcf.gz; do
    sample=$(basename $vcf .vcf.gz)
    echo "Processando $sample..."
    
    python3 tools/vcf_to_frog.py \
        "$vcf" \
        tools/aisnps_55_list.txt \
        "input/${sample}.txt"
done

# Combinar todos em um arquivo
echo "Individual|$(head -1 tools/aisnps_55_list.txt | tr '\n' '|')" > input/all_samples.txt
for file in input/*.txt; do
    tail -n +2 "$file" >> input/all_samples.txt
done

# Analisar todas as amostras juntas
./run.sh
```

## Exemplo 5: Pipeline completo para projeto de pesquisa

```bash
#!/bin/bash
# pipeline_ancestry.sh - Pipeline completo de análise de ancestralidade

set -e

# 1. Baixar dados de referência do 1000 Genomes (primeira vez apenas)
if [ ! -f "input/1000genomes_reference.txt" ]; then
    echo "Baixando dados do 1000 Genomes..."
    ./tools/extract_snps_from_1000genomes.sh \
        -o input/1000genomes_reference.txt \
        -k
fi

# 2. Processar amostras do estudo
echo "Processando amostras do estudo..."
for vcf in study_data/*.vcf.gz; do
    sample=$(basename $vcf .vcf.gz)
    
    ./tools/extract_snps_from_wgs.sh \
        -i "$vcf" \
        -t vcf \
        -o "input/study_${sample}.txt"
done

# 3. Executar análise para cada amostra
for sample_file in input/study_*.txt; do
    sample=$(basename $sample_file .txt)
    
    # Configurar
    cat > FROGAncestryCalc.properties << CONFIG
homePath=.
inputFilename=$(basename $sample_file)
panelInfo=55AI
CONFIG
    
    # Executar
    echo "Analisando $sample..."
    ./run.sh
    
    # Salvar resultados
    mkdir -p results/$sample
    cp output/*_likelihood.txt results/$sample/
    cp output/*_rankOrder.txt results/$sample/
    cp output/*_orderOfMag.txt results/$sample/
done

echo "Pipeline completo!"
```

## Exemplo 6: Verificar qualidade dos dados extraídos

```bash
# Após conversão, verificar quantos SNPs foram capturados
input_file="input/my_sample.txt"

# Contar SNPs no cabeçalho
n_snps=$(head -1 "$input_file" | tr '|' '\n' | tail -n +2 | wc -l)
echo "SNPs encontrados: $n_snps / 55"

# Verificar SNPs faltantes
head -1 "$input_file" | tr '|' '\n' | tail -n +2 > found_snps.txt
comm -23 <(sort tools/aisnps_55_list.txt) <(sort found_snps.txt) > missing_snps.txt

echo "SNPs faltantes:"
cat missing_snps.txt

# Contar genótipos NN (faltantes) por amostra
tail -n +2 "$input_file" | while IFS='|' read -r sample rest; do
    nn_count=$(echo "$rest" | tr '|' '\n' | grep -c "NN" || true)
    echo "$sample: $nn_count genótipos faltantes"
done
```

## Exemplo 7: Filtrar amostras por qualidade antes da conversão

```bash
# Aplicar filtros de qualidade ao VCF antes de converter
input_vcf="raw_data.vcf.gz"
filtered_vcf="filtered_data.vcf.gz"

# Filtrar: QUAL >= 30, DP >= 10, GQ >= 20
bcftools view -i 'QUAL>=30' "$input_vcf" | \
    bcftools view -i 'FORMAT/DP>=10' | \
    bcftools view -i 'FORMAT/GQ>=20' -Oz -o "$filtered_vcf"

bcftools index "$filtered_vcf"

# Agora converter
python3 tools/vcf_to_frog.py \
    "$filtered_vcf" \
    tools/aisnps_55_list.txt \
    input/high_quality_samples.txt
```

## Recursos Adicionais

### Obter genoma de referência

```bash
# GRCh38/hg38 (versão recomendada - usada pelo 1000 Genomes High Coverage)
wget http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# GRCh37/hg19 (versão antiga)
wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
```

### Converter entre builds do genoma

```bash
# Se seus dados estão em hg19 mas precisa hg38
# Use UCSC liftOver ou CrossMap

# Com CrossMap
pip install crossmap
CrossMap.py vcf hg19ToHg38.chain.gz input_hg19.vcf.gz hg38.fa output_hg38.vcf
```

### Verificar build do VCF

```bash
# Checar contig no cabeçalho
bcftools view -h sample.vcf.gz | grep "^##contig"

# GRCh37: contigs geralmente não têm "chr" prefix (1, 2, 3...)
# GRCh38: pode ter "chr" prefix (chr1, chr2, chr3...)
```

## Troubleshooting Comum

### Problema: "No SNPs found"

```bash
# Verificar IDs no VCF
bcftools query -f '%ID\n' sample.vcf.gz | head

# Se não tiver rsIDs, anotar com dbSNP
bcftools annotate -a dbSNP.vcf.gz -c ID sample.vcf.gz -Oz -o annotated.vcf.gz
```

### Problema: "bcftools not found"

```bash
# Instalar via conda
conda install -c bioconda bcftools samtools

# Ou via apt (Ubuntu/Debian)
sudo apt-get install bcftools samtools
```

### Problema: Muitos genótipos NN

- Verifique cobertura do sequenciamento
- Considere imputation de genótipos
- Use filtros de qualidade mais brandos

Para mais informações, consulte `tools/README.md`
