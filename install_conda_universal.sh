#!/bin/bash
# install_conda_universal.sh - Instala conda/mamba em qualquer máquina Linux

set -euo pipefail

echo "🐍 Instalador Universal do Conda/Mamba"
echo "======================================="

# Detecta arquitetura
ARCH=$(uname -m)
OS=$(uname -s)

echo "🖥️  Sistema detectado: $OS $ARCH"

# Define diretório de instalação baseado no usuário
USER_HOME="$HOME"
CONDA_DIR="$USER_HOME/miniforge3"

echo "📁 Diretório de instalação: $CONDA_DIR"

# Verifica se já está instalado
if [ -d "$CONDA_DIR" ]; then
    echo "✅ Conda já instalado em: $CONDA_DIR"
    echo "🔄 Para reinstalar, remova o diretório: rm -rf $CONDA_DIR"
    exit 0
fi

# URLs do Miniforge baseadas na arquitetura
case "$ARCH" in
    x86_64)
        if [ "$OS" = "Linux" ]; then
            MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
        else
            echo "❌ Sistema não suportado: $OS $ARCH"
            exit 1
        fi
        ;;
    aarch64|arm64)
        if [ "$OS" = "Linux" ]; then
            MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"
        else
            echo "❌ Sistema não suportado: $OS $ARCH"
            exit 1
        fi
        ;;
    *)
        echo "❌ Arquitetura não suportada: $ARCH"
        echo "💡 Visite: https://github.com/conda-forge/miniforge#download"
        exit 1
        ;;
esac

echo "📥 Baixando Miniforge..."
echo "🔗 URL: $MINIFORGE_URL"

# Baixa instalador
INSTALLER="/tmp/miniforge_installer.sh"
curl -fsSL "$MINIFORGE_URL" -o "$INSTALLER"

echo "✅ Download concluído"

# Torna executável e instala
chmod +x "$INSTALLER"

echo "🔧 Instalando Miniforge..."
echo "📁 Destino: $CONDA_DIR"

# Instala em modo silencioso
bash "$INSTALLER" -b -p "$CONDA_DIR"

echo "✅ Miniforge instalado!"

# Limpa instalador
rm -f "$INSTALLER"

# Inicializa conda
echo "🔧 Inicializando conda..."
"$CONDA_DIR/bin/conda" init bash

echo "🐍 Instalando mamba..."
"$CONDA_DIR/bin/conda" install -n base -c conda-forge mamba -y

echo ""
echo "✅ Instalação concluída!"
echo "📁 Conda instalado em: $CONDA_DIR"
echo ""
echo "🚀 Próximos passos:"
echo "   1. Reinicie o terminal ou execute: source ~/.bashrc"
echo "   2. Execute: ./install_genomics_env.sh"
echo "   3. Execute: source start_genomics.sh"
echo ""
echo "💡 Para testar:"
echo "   conda --version"
echo "   mamba --version"
