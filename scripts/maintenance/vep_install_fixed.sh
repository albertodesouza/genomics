#!/bin/bash
# vep_install_fixed.sh - Instalador VEP corrigido e simplificado

echo "🧬 Instalador VEP Corrigido - Versão Estável"
echo "==========================================="

# Configurações
SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
VEPCACHE="/dados/vep_cache"
VEP_DIR="$HOME/ensembl-vep"

echo "⚙️  Configurações:"
echo "   Espécie: $SPECIES"
echo "   Assembly: $ASSEMBLY" 
echo "   Cache: $VEPCACHE"
echo "   Diretório: $VEP_DIR"

# Verifica ambiente conda
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "❌ Ambiente conda não ativo"
    echo "💡 Execute: source start_genomics_universal.sh"
    exit 1
fi

echo "✅ Ambiente: $CONDA_DEFAULT_ENV"

# Cria diretórios
mkdir -p "$VEPCACHE"

echo ""
echo "📦 1. Instalando dependências Perl..."

# Instala dependências essenciais via conda
conda install -c bioconda -c conda-forge -y \
    ensembl-vep \
    perl-list-moreutils \
    perl-set-intervaltree \
    perl-dbi \
    perl-archive-zip || {
    
    echo "⚠️  Instalação via conda falhou, tentando método manual..."
}

echo ""
echo "🧪 2. Testando VEP básico..."

if command -v vep >/dev/null 2>&1; then
    echo "✅ VEP encontrado via conda: $(which vep)"
    VEP_VERSION=$(vep 2>&1 | head -10 | grep -E "ensembl-vep.*[0-9]" || echo "versão conda")
    echo "📋 Versão: $VEP_VERSION"
else
    echo "❌ VEP não encontrado via conda"
    echo "🔄 Tentando instalação manual do GitHub..."
    
    # Remove diretório anterior se existir
    if [ -d "$VEP_DIR" ]; then
        rm -rf "$VEP_DIR"
    fi
    
    # Clona versão estável
    echo "📥 Clonando VEP release/114..."
    git clone -b release/114 --depth 1 https://github.com/Ensembl/ensembl-vep.git "$VEP_DIR"
    
    cd "$VEP_DIR"
    
    # Instalação mínima
    echo "🔧 Instalação mínima do VEP..."
    chmod +x INSTALL.pl
    ./INSTALL.pl --AUTO a --NO_BIOPERL --NO_HTSLIB --NO_TEST --NO_UPDATE
    
    # Configura PATH
    export PATH="$VEP_DIR:$PATH"
    export PERL5LIB="$VEP_DIR:$PERL5LIB"
    
    # Adiciona ao ambiente conda
    CONDA_PREFIX="${CONDA_PREFIX:-$HOME/miniforge3/envs/genomics}"
    if [ -d "$CONDA_PREFIX/etc/conda/activate.d" ]; then
        cat > "$CONDA_PREFIX/etc/conda/activate.d/vep_manual.sh" << EOF
export PATH="$VEP_DIR:\$PATH"
export PERL5LIB="$VEP_DIR:\$PERL5LIB"
EOF
        echo "✅ PATH configurado para conda"
    fi
fi

echo ""
echo "💾 3. Instalando cache..."

# Tenta instalar cache
if command -v vep >/dev/null 2>&1; then
    echo "🔧 Baixando cache para $SPECIES $ASSEMBLY..."
    
    # Método 1: via vep_install se disponível
    if command -v vep_install >/dev/null 2>&1; then
        echo "📦 Usando vep_install (método preferido)..."
        vep_install -a cf -s "$SPECIES" -y "$ASSEMBLY" -c "$VEPCACHE" --NO_BIOPERL || {
            echo "⚠️  vep_install falhou, tentando INSTALL.pl..."
        }
    fi
    
    # Método 2: via INSTALL.pl
    if [ -f "$VEP_DIR/INSTALL.pl" ]; then
        echo "📦 Usando INSTALL.pl..."
        cd "$VEP_DIR"
        ./INSTALL.pl --AUTO cf --SPECIES "$SPECIES" --ASSEMBLY "$ASSEMBLY" --CACHEDIR "$VEPCACHE" --NO_BIOPERL || {
            echo "⚠️  Cache automático falhou"
        }
    fi
else
    echo "❌ VEP não encontrado para instalar cache"
fi

echo ""
echo "🧪 4. Teste final..."

# Teste completo
if command -v vep >/dev/null 2>&1; then
    echo "✅ VEP disponível: $(which vep)"
    
    # Testa help
    if vep --help >/dev/null 2>&1; then
        echo "✅ VEP funcional"
        
        # Verifica cache
        if [ -d "$VEPCACHE/$SPECIES" ]; then
            CACHE_FOUND=$(find "$VEPCACHE/$SPECIES" -name "*$ASSEMBLY*" -type d | head -1)
            if [ -n "$CACHE_FOUND" ]; then
                echo "✅ Cache encontrado: $CACHE_FOUND"
            else
                echo "⚠️  Cache não encontrado"
            fi
        fi
        
        echo ""
        echo "🎉 Instalação VEP concluída!"
        echo "📋 Para usar:"
        echo "   vep --cache --offline --species $SPECIES --assembly $ASSEMBLY --dir_cache $VEPCACHE"
        
    else
        echo "❌ VEP com problemas de dependências"
        echo "💡 Tente:"
        echo "   conda install -c bioconda perl-list-moreutils"
        echo "   source start_genomics_universal.sh"
    fi
else
    echo "❌ VEP não instalado"
    echo "💡 Tente instalação manual:"
    echo "   conda install -c bioconda ensembl-vep"
fi

echo ""
echo "💡 Se houver problemas, use o instalador original:"
echo "   source vep_install.sh"
