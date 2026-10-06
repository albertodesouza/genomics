#!/bin/bash
# vep_install_latest.sh - Instala a última versão do VEP diretamente do GitHub Ensembl

set -euo pipefail

echo "🧬 Instalador VEP - Última Versão do GitHub Ensembl"
echo "================================================="

# Configurações
SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
VEPCACHE="/dados/vep_cache"
VEP_DIR="$HOME/ensembl-vep"
CONDA_ENV="genomics"

echo "⚙️  Configurações:"
echo "   Espécie: $SPECIES"
echo "   Assembly: $ASSEMBLY"
echo "   Cache: $VEPCACHE"
echo "   Diretório VEP: $VEP_DIR"
echo ""

# Verifica se conda está ativo
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]; then
    echo "❌ Ambiente conda '$CONDA_ENV' não está ativo"
    echo "💡 Execute: source start_genomics_universal.sh"
    exit 1
fi

echo "✅ Ambiente conda ativo: $CONDA_DEFAULT_ENV"

# Cria diretórios
mkdir -p "$VEPCACHE"
mkdir -p "$(dirname "$VEP_DIR")"

echo ""
echo "📥 1. Baixando última versão do VEP do GitHub..."

# Remove instalação anterior se existir
if [ -d "$VEP_DIR" ]; then
    echo "🧹 Removendo instalação anterior..."
    rm -rf "$VEP_DIR"
fi

# Clona repositório oficial do Ensembl
echo "📦 Clonando repositório oficial do Ensembl VEP..."
git clone https://github.com/Ensembl/ensembl-vep.git "$VEP_DIR"

cd "$VEP_DIR"

# Verifica versões disponíveis
echo ""
echo "🔍 Verificando versões disponíveis..."
git tag --sort=-version:refname | head -10

# Pega a última release (não beta/rc)
LATEST_VERSION=$(git tag --sort=-version:refname | grep -E '^[0-9]+\.[0-9]+(\.[0-9]+)?$' | head -1)
if [ -z "$LATEST_VERSION" ]; then
    echo "⚠️  Não foi possível detectar versão estável, usando release/114"
    LATEST_VERSION="release/114"
fi

echo "🎯 Usando versão: $LATEST_VERSION"

# Checkout da versão mais recente
git checkout "$LATEST_VERSION"

echo ""
echo "🔧 2. Instalando VEP..."

# Verifica dependências Perl
echo "🐪 Verificando dependências Perl..."
perl -MDBI -e 1 || {
    echo "❌ Perl DBI ausente - instalando..."
    conda install -c bioconda perl-dbi -y
}

perl -MArchive::Zip -e 1 2>/dev/null || {
    echo "📦 Instalando Archive::Zip..."
    conda install -c conda-forge perl-archive-zip -y
}

# Instala VEP
echo "⚡ Executando instalação do VEP..."
echo "💡 Isso pode demorar 10-20 minutos..."

# Torna executável
chmod +x INSTALL.pl

# Executa instalação com parâmetros otimizados
./INSTALL.pl \
    --AUTO a \
    --SPECIES "$SPECIES" \
    --ASSEMBLY "$ASSEMBLY" \
    --CACHEDIR "$VEPCACHE" \
    --NO_BIOPERL \
    --NO_HTSLIB \
    --NO_TEST \
    --VERBOSE

echo ""
echo "🔗 3. Configurando PATH..."

# Adiciona VEP ao PATH do ambiente conda
VEP_BIN="$VEP_DIR"
CONDA_PREFIX="${CONDA_PREFIX:-$HOME/miniforge3/envs/$CONDA_ENV}"

if [ -d "$CONDA_PREFIX" ]; then
    # Cria script de ativação
    mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
    cat > "$CONDA_PREFIX/etc/conda/activate.d/vep_path.sh" << EOF
#!/bin/bash
# Adiciona VEP ao PATH
export PATH="$VEP_BIN:\$PATH"
export PERL5LIB="$VEP_DIR:\$PERL5LIB"
EOF
    
    # Ativa imediatamente
    export PATH="$VEP_BIN:$PATH"
    export PERL5LIB="$VEP_DIR:$PERL5LIB"
    
    echo "✅ PATH configurado para ambiente conda"
else
    echo "⚠️  Configuração manual necessária:"
    echo "   export PATH=\"$VEP_BIN:\$PATH\""
    echo "   export PERL5LIB=\"$VEP_DIR:\$PERL5LIB\""
fi

echo ""
echo "🧪 4. Testando instalação..."

# Testa VEP
if command -v vep >/dev/null 2>&1; then
    VEP_VERSION=$(vep --help 2>&1 | head -20 | grep -E "(ensembl-vep|version)" | head -1 || echo "instalado")
    echo "✅ VEP instalado com sucesso!"
    echo "📋 Versão: $VEP_VERSION"
    echo "📁 Localização: $(which vep)"
else
    echo "❌ VEP não encontrado no PATH"
    echo "💡 Tente reiniciar o terminal ou executar:"
    echo "   source start_genomics_universal.sh"
fi

echo ""
echo "📊 5. Verificando cache..."

# Verifica cache baixado
if [ -d "$VEPCACHE/$SPECIES" ]; then
    echo "✅ Cache encontrado:"
    ls -la "$VEPCACHE/$SPECIES"/ | head -5
    
    # Procura diretório da versão
    CACHE_DIRS=$(find "$VEPCACHE/$SPECIES" -maxdepth 1 -type d -name "*$ASSEMBLY*" 2>/dev/null || true)
    if [ -n "$CACHE_DIRS" ]; then
        echo "✅ Cache para $ASSEMBLY encontrado:"
        echo "$CACHE_DIRS"
    else
        echo "⚠️  Cache para $ASSEMBLY não encontrado"
        echo "💡 Execute manualmente:"
        echo "   vep_install -a cf -s $SPECIES -y $ASSEMBLY -c $VEPCACHE --NO_BIOPERL"
    fi
else
    echo "❌ Cache não encontrado em $VEPCACHE"
fi

echo ""
echo "🎯 Instalação Concluída!"
echo "======================"
echo "📋 Resumo:"
echo "   • VEP: $(which vep 2>/dev/null || echo 'não encontrado no PATH')"
echo "   • Versão: $LATEST_VERSION"
echo "   • Cache: $VEPCACHE"
echo "   • Espécie: $SPECIES"
echo "   • Assembly: $ASSEMBLY"
echo ""
echo "💡 Para usar:"
echo "   vep --help"
echo "   vep --cache --offline --species $SPECIES --assembly $ASSEMBLY --dir_cache $VEPCACHE"
echo ""
echo "🔄 Se houver problemas, reinicie o terminal e execute:"
echo "   source start_genomics_universal.sh"
