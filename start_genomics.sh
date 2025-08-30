#!/usr/bin/env bash
# Ativa o conda de forma local e entra no env "genomics"
# Uso recomendado:  source ~/bin/start_genomics.sh

# set -euo pipefail

CONDA_BASE="~/miniforge3"

# Inicializa o conda sem precisar do bloco no .bashrc
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1090
  . "$CONDA_BASE/etc/profile.d/conda.sh"
else
  export PATH="$CONDA_BASE/bin:$PATH"
fi

conda activate genomics

