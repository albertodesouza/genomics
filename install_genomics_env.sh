#!/usr/bin/env bash
# Verbose installer for the "genomics" environment + optional VEP cache download
set -euo pipefail

#######################################
# Defaults
#######################################
ENV_NAME="genomics"
INSTALL_VEP_CACHE=0
VEP_SPECIES="homo_sapiens"
VEP_ASSEMBLY="GRCh38"
VEP_CACHE_DIR="${HOME}/.vep"
VEP_FORK="8"
DEBUG=0

PACKAGES=(
  python=3.11 pyyaml rich tqdm humanfriendly psutil
  sra-tools fastqc multiqc cutadapt
  bwa bwa-mem2 minimap2 samtools picard gatk4 bedtools
  gffread ensembl-vep bcftools
  hisat2 stringtie gffcompare wget
  mosdepth seqtk
)

#######################################
# Helpers
#######################################
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  --env-name NAME            Conda env name (default: ${ENV_NAME})
  --install-vep-cache        Download and install the VEP cache
  --vep-species SPECIES      VEP species (default: ${VEP_SPECIES})
  --vep-assembly ASSEMBLY    VEP assembly (GRCh38|GRCh37; default: ${VEP_ASSEMBLY})
  --vep-cache-dir DIR        VEP cache directory (default: ${VEP_CACHE_DIR})
  --vep-fork N               Parallel download/extract for VEP cache (default: ${VEP_FORK})
  --debug                    Enable 'set -x' for command tracing
  -h, --help                 Show this help
EOF
}

norm_assembly() {
  local a="${1^^}"
  if [[ "$a" == *"GRCH38"* ]]; then echo "GRCh38"; return; fi
  if [[ "$a" == *"GRCH37"* || "$a" == *"HG19"* ]]; then echo "GRCh37"; return; fi
  echo "GRCh38"
}

die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo -e "\n\033[1;34m==>\033[0m $*"; }

#######################################
# Arg parsing
#######################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-name)          ENV_NAME="${2:?}"; shift 2 ;;
    --install-vep-cache) INSTALL_VEP_CACHE=1; shift ;;
    --vep-species)       VEP_SPECIES="${2:?}"; shift 2 ;;
    --vep-assembly)      VEP_ASSEMBLY="$(norm_assembly "${2:?}")"; shift 2 ;;
    --vep-cache-dir)     VEP_CACHE_DIR="${2:?}"; shift 2 ;;
    --vep-fork)          VEP_FORK="${2:?}"; shift 2 ;;
    --debug)             DEBUG=1; shift ;;
    -h|--help)           usage; exit 0 ;;
    *) die "Unknown option: $1 (use --help)";;
  esac
done

[[ $DEBUG -eq 1 ]] && set -x

#######################################
# Pre-checks (no /dev/null)
#######################################
info "Checking for conda..."
command -v conda || die "conda not found in PATH."

info "Checking for mamba (optional speed-up)..."
if command -v mamba; then
  echo "mamba found: $(command -v mamba)"
else
  info "Installing mamba in base env (full output below)"
  conda install -n base -c conda-forge -y mamba
fi

# Decide package manager
if command -v mamba; then
  PM="mamba"
else
  PM="conda"
fi
echo "Package manager: ${PM}"

#######################################
# Create or update environment (verbose)
#######################################
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  info "Environment '${ENV_NAME}' exists. Ensuring required packages..."
  "$PM" install -y -n "${ENV_NAME}" -c bioconda -c conda-forge --strict-channel-priority "${PACKAGES[@]}"
else
  info "Creating environment '${ENV_NAME}'..."
  "$PM" create  -y -n "${ENV_NAME}" -c bioconda -c conda-forge --strict-channel-priority "${PACKAGES[@]}"
  echo "To activate later: conda activate ${ENV_NAME}"
fi

#######################################
# VEP cache (optional) — fully verbose
#######################################
if [[ $INSTALL_VEP_CACHE -eq 1 ]]; then
  info "Preparing VEP cache → species='${VEP_SPECIES}', assembly='${VEP_ASSEMBLY}', dir='${VEP_CACHE_DIR}', fork=${VEP_FORK}"

  mkdir -p "${VEP_CACHE_DIR}"
  [[ -w "${VEP_CACHE_DIR}" ]] || die "Cache dir '${VEP_CACHE_DIR}' is not writable."

  # Runner that does NOT capture output
  if command -v mamba; then
    RUNNER=(mamba run -n "${ENV_NAME}")
  else
    RUNNER=(conda run --no-capture-output -n "${ENV_NAME}")
  fi

  info "Checking for 'vep_install' inside environment..."
  if "${RUNNER[@]}" bash -lc 'command -v vep_install'; then
    info "Using 'vep_install' (modern installer)."
    "${RUNNER[@]}" vep_install \
      -a cf \
      -s "${VEP_SPECIES}" \
      -y "${VEP_ASSEMBLY}" \
      -c "${VEP_CACHE_DIR}" \
      --AUTO c \
      --NO_UPDATE \
      --CONVERT \
      --fork "${VEP_FORK}"
  else
    info "'vep_install' not found; searching for legacy INSTALL.pl under \$CONDA_PREFIX/share ..."
    # This will print errors if the glob doesn't match; that's intentional (no suppression).
    INSTALL_PL_PATH="$("${RUNNER[@]}" bash -lc 'ls -d "${CONDA_PREFIX}"/share/ensembl-vep-*/INSTALL.pl | head -n1 || true')"
    echo "INSTALL.pl path: ${INSTALL_PL_PATH:-<not found>}"
    [[ -n "${INSTALL_PL_PATH}" ]] || die "Could not locate INSTALL.pl under \$CONDA_PREFIX/share."

    info "Running legacy INSTALL.pl (full output below)."
    "${RUNNER[@]}" perl "${INSTALL_PL_PATH}" \
      -a cf \
      -s "${VEP_SPECIES}" \
      -y "${VEP_ASSEMBLY}" \
      -c "${VEP_CACHE_DIR}" \
      --AUTO c \
      --NO_UPDATE \
      --CONVERT \
      --fork "${VEP_FORK}"
  fi

  info "VEP cache finished. Showing directory tree:"
  # Show up to two levels; if 'tree' missing, fall back to ls -lR
  if command -v tree; then
    tree -L 2 "${VEP_CACHE_DIR}" || true
  else
    ls -lR "${VEP_CACHE_DIR}" || true
  fi
  echo "Point your YAML 'params.vep_dir_cache' to: ${VEP_CACHE_DIR}"
fi

info "All done."

