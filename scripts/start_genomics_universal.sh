#!/usr/bin/env bash
# start_genomics_universal.sh - Ativa ambiente genomics em qualquer máquina
# Uso: source start_genomics_universal.sh

# Detecta automaticamente onde o conda está instalado
detect_conda_base() {
    # Lista de locais comuns para conda/miniforge
    local conda_paths=(
        "$HOME/miniforge3"
        "$HOME/miniconda3" 
        "$HOME/anaconda3"
        "/opt/conda"
        "/opt/miniforge3"
        "/usr/local/conda"
        "/usr/local/miniforge3"
        "/home/lume2/miniforge3"  # seu path específico
    )
    
    # Verifica cada local
    for path in "${conda_paths[@]}"; do
        if [ -d "$path" ] && [ -f "$path/etc/profile.d/conda.sh" ]; then
            echo "$path"
            return 0
        fi
    done
    
    # Se não encontrou, tenta detectar via comando conda
    if command -v conda >/dev/null 2>&1; then
        local conda_exe=$(which conda)
        local conda_base=$(dirname $(dirname "$conda_exe"))
        if [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
            echo "$conda_base"
            return 0
        fi
    fi
    
    return 1
}

echo "🐍 Inicializando ambiente genomics..."

# Detecta conda automaticamente
CONDA_BASE=$(detect_conda_base)

if [ -z "$CONDA_BASE" ]; then
    echo "❌ Conda não encontrado!"
    echo ""
    echo "💡 Para instalar:"
    echo "   ./install_conda_universal.sh"
    echo ""
    echo "🔍 Locais verificados:"
    echo "   $HOME/miniforge3"
    echo "   $HOME/miniconda3"
    echo "   $HOME/anaconda3"
    echo "   /opt/conda"
    echo "   /opt/miniforge3"
    echo "   /home/lume2/miniforge3"
    return 1
fi

echo "✅ Conda encontrado: $CONDA_BASE"

# Inicializa conda
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1090
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    export PATH="$CONDA_BASE/bin:$PATH"
fi

# Verifica se ambiente genomics existe
if ! conda env list | grep -q "^genomics "; then
    echo "❌ Ambiente 'genomics' não encontrado!"
    echo ""
    echo "💡 Para criar:"
    echo "   ./install_genomics_env.sh"
    echo ""
    echo "🔍 Ambientes disponíveis:"
    conda env list
    return 1
fi

# Ativa ambiente
echo "🔧 Ativando ambiente genomics..."
conda activate genomics

# Evita misturar bibliotecas CUDA do sistema com as bibliotecas empacotadas
# pelo PyTorch dentro do ambiente conda. Isso quebra a inicialização CUDA em
# alguns setups quando /usr/local/cuda/lib64 entra no LD_LIBRARY_PATH.
if [ -n "$LD_LIBRARY_PATH" ]; then
    CLEAN_LD_LIBRARY_PATH=$(python3 - <<'PY'
import os

parts = os.environ.get("LD_LIBRARY_PATH", "").split(":")
filtered = [
    p for p in parts
    if p and not p.startswith("/usr/local/cuda")
]
print(":".join(filtered))
PY
)
    export LD_LIBRARY_PATH="$CLEAN_LD_LIBRARY_PATH"
fi
unset CUDA_LIBS

# Verifica se ativação funcionou
if [ "$CONDA_DEFAULT_ENV" = "genomics" ]; then
    echo "✅ Ambiente genomics ativo!"
    echo "🧬 Ferramentas disponíveis:"
    
    # Verifica ferramentas principais
    tools=("bwa" "bwa-mem2" "samtools" "bcftools" "vep" "seqtk")
    for tool in "${tools[@]}"; do
        if command -v "$tool" >/dev/null 2>&1; then
            version=$($tool --version 2>/dev/null | head -1 || $tool 2>&1 | head -1 || echo "instalado")
            echo "   ✅ $tool: $version"
        else
            echo "   ❌ $tool: não encontrado"
        fi
    done
    
    echo ""
    echo "🎯 Sistema pronto para análise genômica!"
    echo "📁 Diretório atual: $(pwd)"
else
    echo "❌ Falha ao ativar ambiente genomics"
    return 1
fi
