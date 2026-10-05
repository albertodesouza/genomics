#!/usr/bin/env bash
# install.sh - one-command install of the genomics workspace for the visualizer (and optionally
# training and a local AlphaGenome server). Safe to re-run: existing pieces are kept.
#
#   scripts/env/install.sh                         # conda env "genomics": visualizer + import + AlphaGenome client
#   scripts/env/install.sh --training              # + PyTorch/scikit-learn for training and the Perturbation Lab
#   scripts/env/install.sh --alphagenome-server    # + AlphaGenome on this machine's NVIDIA GPU (separate env)
#   scripts/env/install.sh --no-conda              # plain venv in .venv (install bcftools/samtools yourself)
#
# What each feature needs, with disk and hardware sizes: docs/getting-started/requirements.md
set -euo pipefail

ENV_NAME="genomics"
PREFIX=""
PYTHON_VERSION="3.11"
TRAINING=0
AG_SERVER=0
USE_CONDA=1
TORCH_INDEX="${TORCH_INDEX:-}"

usage() { sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV_NAME="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --python) PYTHON_VERSION="$2"; shift 2 ;;
    --training) TRAINING=1; shift ;;
    --alphagenome-server) AG_SERVER=1; shift ;;
    --no-conda) USE_CONDA=0; shift ;;
    -h|--help) usage 0 ;;
    *) echo "unknown option: $1" >&2; usage 2 ;;
  esac
done

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ARCH="$(uname -m)"
step() { printf '\n==> %s\n' "$*"; }

EXTRAS="visualizer"
[[ "${TRAINING}" == 1 ]] && EXTRAS="visualizer,genotype"

if [[ "${USE_CONDA}" == 1 ]]; then
  CONDA="${CONDA_EXE:-$(command -v mamba || command -v conda || true)}"
  if [[ -z "${CONDA}" ]]; then
    for base in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda /opt/miniforge3; do
      [[ -x "${base}/bin/conda" ]] && CONDA="${base}/bin/conda" && break
    done
  fi
  if [[ -z "${CONDA}" ]]; then
    echo "conda not found. Install Miniforge first (scripts/env/install_conda_universal.sh), or re-run with --no-conda." >&2
    exit 1
  fi
  if [[ -n "${PREFIX}" ]]; then
    TARGET=(-p "${PREFIX}"); ENV_DIR="${PREFIX}"
  else
    TARGET=(-n "${ENV_NAME}"); ENV_DIR="$("${CONDA}" info --base)/envs/${ENV_NAME}"
  fi
  PY="${ENV_DIR}/bin/python"
  # bcftools/samtools/htslib: VCF import, new prediction windows and the training-axis alignment.
  TOOLS=(bcftools samtools htslib)
  if [[ ! -x "${PY}" ]]; then
    step "Creating conda env ${ENV_DIR} (Python ${PYTHON_VERSION}, bcftools, samtools)"
    "${CONDA}" create -y "${TARGET[@]}" -c conda-forge -c bioconda "python=${PYTHON_VERSION}" pip "${TOOLS[@]}"
  elif [[ ! -x "${ENV_DIR}/bin/bcftools" || ! -x "${ENV_DIR}/bin/samtools" ]]; then
    step "Adding bcftools and samtools to ${ENV_DIR}"
    "${CONDA}" install -y "${TARGET[@]}" -c conda-forge -c bioconda "${TOOLS[@]}"
  else
    step "Using existing conda env ${ENV_DIR}"
  fi
else
  ENV_DIR="${PREFIX:-${REPO}/.venv}"
  PY="${ENV_DIR}/bin/python"
  if [[ ! -x "${PY}" ]]; then
    step "Creating venv ${ENV_DIR}"
    python3 -m venv "${ENV_DIR}"
  fi
  for tool in bcftools samtools; do
    command -v "${tool}" >/dev/null 2>&1 || echo "note: ${tool} is not on PATH; dataset import needs it (apt install bcftools samtools, or conda)." >&2
  done
fi

step "Installing genomics with extras [${EXTRAS}]"
"${PY}" -m pip install --upgrade pip
if [[ "${TRAINING}" == 1 && -z "${TORCH_INDEX}" && "${ARCH}" == "aarch64" ]] && command -v nvidia-smi >/dev/null 2>&1; then
  # PyPI's aarch64 torch wheels are CPU-only; NVIDIA ARM machines (DGX Spark, GH200) need PyTorch's CUDA index.
  TORCH_INDEX="https://download.pytorch.org/whl/cu128"
fi
if [[ -n "${TORCH_INDEX}" ]]; then
  "${PY}" -m pip install torch --index-url "${TORCH_INDEX}"
fi
"${PY}" -m pip install -e "${REPO}[${EXTRAS}]"

if [[ "${AG_SERVER}" == 1 ]]; then
  step "Setting up the local AlphaGenome server (separate environment, CUDA jax)"
  "${ENV_DIR}/bin/genomics" alphagenome server setup || echo "AlphaGenome server setup is incomplete; see the messages above." >&2
fi

step "Checking what this machine can run"
"${ENV_DIR}/bin/genomics" doctor || true

echo
if [[ "${USE_CONDA}" == 1 && -z "${PREFIX}" ]]; then
  echo "Done. Activate the environment, then start the visualizer:"
  echo "  conda activate ${ENV_NAME}        # or: source scripts/env/start_genomics_universal.sh"
else
  echo "Done. Activate the environment, then start the visualizer:"
  echo "  source ${ENV_DIR}/bin/activate     # conda: conda activate ${ENV_DIR}"
fi
echo "  genomics visualize --open"
