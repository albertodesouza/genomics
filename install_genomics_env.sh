#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="genomics"

# Mamba sobre conda (se faltar)
if ! command -v mamba >/dev/null 2>&1; then
  conda install -n base -c conda-forge -y mamba
fi

# Se o ambiente não existir, cria; se existir, apenas garante os pacotes
if ! conda env list | grep -q "^${ENV_NAME}\b"; then
  mamba create -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
    python=3.11 pyyaml rich tqdm humanfriendly psutil \
    sra-tools fastqc multiqc cutadapt \
    bwa-mem2 minimap2 samtools picard gatk4 bedtools \
    gffread ensembl-vep bcftools \
    hisat2 stringtie gffcompare wget \
    mosdepth seqtk
  echo "Ative com: conda activate ${ENV_NAME}"
else
  # Já está criado: atualiza/garante pacotes, podendo rodar DENTRO do ambiente
  mamba install -y -c bioconda -c conda-forge \
    pyyaml rich tqdm humanfriendly psutil \
    sra-tools fastqc multiqc cutadapt \
    bwa-mem2 minimap2 samtools picard gatk4 bedtools \
    gffread ensembl-vep bcftools \
    hisat2 stringtie gffcompare wget \
    mosdepth seqtk
fi

