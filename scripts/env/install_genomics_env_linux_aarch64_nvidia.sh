#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${ENV_NAME:-genomics}"
PYTHON_VERSION="${PYTHON_VERSION:-3.10}"
TORCH_VERSION="${TORCH_VERSION:-2.8.0}"
TORCH_CUDA_TAG="${TORCH_CUDA_TAG:-cu128}"

ARCH="$(uname -m)"
OS="$(uname -s)"

if [[ "${OS}" != "Linux" ]]; then
  printf "Este script suporta apenas Linux. Detectado: %s\n" "${OS}" >&2
  exit 1
fi

if [[ "${ARCH}" != "aarch64" ]]; then
  printf "Este script e especifico para linux-aarch64. Detectado: %s\n" "${ARCH}" >&2
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  printf "conda nao encontrado no PATH. Instale Miniconda/Anaconda primeiro.\n" >&2
  exit 1
fi

if ! command -v nvidia-smi >/dev/null 2>&1; then
  printf "nvidia-smi nao encontrado. Verifique driver NVIDIA na DGX Spark.\n" >&2
  exit 1
fi

if ! command -v mamba >/dev/null 2>&1; then
  conda install -n base -c conda-forge -y mamba
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | grep -q "^${ENV_NAME}\\b"; then
  mamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda python="${PYTHON_VERSION}" perl pip
fi

conda activate "${ENV_NAME}"

printf "\n[1/4] Instalando stack base bioinformatica viavel em linux-aarch64...\n"
mamba install -y -c conda-forge -c bioconda --channel-priority flexible \
  bcftools samtools htslib \
  pyyaml rich tqdm humanfriendly psutil \
  numpy scipy pandas scikit-learn pydantic matplotlib seaborn jupyterlab \
  sra-tools fastqc multiqc cutadapt \
  bwa minimap2 bedtools \
  gffread \
  hisat2 stringtie gffcompare wget \
  mosdepth seqtk pigz \
  perl-dbi perl-app-cpanminus \
  perl-archive-zip perl-json perl-try-tiny

printf "\n[2/4] Instalando dependencias Python auxiliares...\n"
python -m pip install --upgrade pip
python -m pip install \
  xgboost \
  datasets

printf "\n[3/4] Instalando PyTorch para linux-aarch64 + NVIDIA...\n"
python -m pip install --no-cache-dir \
  "torch==${TORCH_VERSION}" \
  --index-url "https://download.pytorch.org/whl/${TORCH_CUDA_TAG}"

printf "\n[4/4] Validando instalacao principal...\n"
python - <<'PY'
import shutil
import torch
import torch._dynamo

print('torch_version=', torch.__version__)
print('cuda_available=', torch.cuda.is_available())
print('cuda_version=', torch.version.cuda)
print('device_count=', torch.cuda.device_count())
print('torch_dynamo_ok=', True)
print('bcftools_in_path=', shutil.which('bcftools') is not None)
PY

printf "\nAmbiente '%s' pronto para linux-aarch64 + NVIDIA.\n" "${ENV_NAME}"
printf "Pacotes omitidos propositalmente neste alvo ARM: gatk4, admixture, plink, plink2, bwa-mem2, picard.\n"
printf "Ative com: conda activate %s\n" "${ENV_NAME}"
