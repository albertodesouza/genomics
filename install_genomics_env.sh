#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="genomics"
PYTHON_VERSION="3.10"

# Instalar mamba no base se não existir
if ! command -v mamba >/dev/null 2>&1; then
  conda install -n base -c conda-forge -y mamba
fi

# Criar ambiente com python primeiro, se não existir
if ! conda env list | grep -q "^${ENV_NAME}\\b"; then
  mamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda python=${PYTHON_VERSION} perl
fi

# Ativar o ambiente
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Instalar pacotes principais, incluindo perl-dbi e ensembl-vep
mamba install -y -c conda-forge -c bioconda --channel-priority flexible \
  bcftools samtools htslib pyyaml rich tqdm humanfriendly psutil \
  sra-tools fastqc multiqc cutadapt \
  bwa bwa-mem2 minimap2 samtools picard gatk4 bedtools \
  gffread ensembl-vep bcftools \
  hisat2 stringtie gffcompare wget \
  mosdepth seqtk pigz \
  perl-dbi perl-app-cpanminus

# Instalar módulo Perl DBI atualizado via cpanminus para evitar erro de versão
cpanm DBI

echo "Ambiente '${ENV_NAME}' pronto. Ative com: conda activate ${ENV_NAME}"

