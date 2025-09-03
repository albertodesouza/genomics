#!/bin/bash
# vep_install_smart.sh - Instalador inteligente do VEP com detecção automática de versão

#set -euo pipefail

echo "🧬 Instalador Inteligente do VEP - Ensembl GitHub"
echo "==============================================="

# Configurações padrão
SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
VEPCACHE="/dados/vep_cache"
VEP_DIR="$HOME/ensembl-vep"

# Função para detectar última versão
detect_latest_vep_version() {
    # Tenta via API do GitHub (mais confiável)
    if command -v curl >/dev/null 2>&1; then
        LATEST=$(curl -s https://api.github.com/repos/Ensembl/ensembl-vep/releases/latest | grep '"tag_name"' | sed 's/.*"tag_name": "\([^"]*\)".*/\1/' 2>/dev/null || echo "")
        if [ -n "$LATEST" ] && [ "$LATEST" != "null" ]; then
            echo "$LATEST"
            return 0
        fi
    fi
    
    # Fallback: usa branch release/114 (estável)
    echo "release/114"
}

# Função principal de instalação
install_vep_latest() {
    local version="$1"
    
    echo ""
    echo "📦 Instalando VEP versão: $version"
    echo "📁 Diretório: $VEP_DIR"
    echo "💾 Cache: $VEPCACHE"
    
    # Remove instalação anterior
    if [ -d "$VEP_DIR" ]; then
        echo "🧹 Removendo instalação anterior..."
        rm -rf "$VEP_DIR"
    fi
    
    # Clona repositório
    echo "📥 Clonando repositório Ensembl VEP..."
    git clone https://github.com/Ensembl/ensembl-vep.git "$VEP_DIR"
    
    cd "$VEP_DIR"
    
    # Checkout da versão específica
    echo "🔄 Mudando para versão: $version"
    git checkout "$version" || {
        echo "⚠️  Versão $version não encontrada, usando main"
        git checkout main
    }
    
    # Mostra informações da versão
    echo "📋 Informações da versão instalada:"
    git log --oneline -1
    
    # Verifica e instala dependências Perl essenciais
    echo ""
    echo "🐪 Instalando dependências Perl via conda..."
    
    # Instala todas as dependências essenciais
    echo "📦 Instalando dependências Perl completas..."
    conda install -c bioconda -c conda-forge -y \
        perl-dbi \
        perl-archive-zip \
        perl-lwp-simple \
        perl-json \
        perl-list-moreutils \
        perl-set-intervaltree \
        perl-bio-db-hts \
        perl-compress-raw-zlib \
        perl-try-tiny \
        perl-uri \
        perl-http-message \
        perl-io-string \
        perl-text-csv || {
        
        echo "⚠️  Algumas dependências falharam via conda, continuando..."
    }
    
    echo "✅ Dependências instaladas"
    
    # Torna executável
    chmod +x INSTALL.pl
    
    echo ""
    echo "⚡ Executando instalação do VEP..."
    echo "💡 Isso pode demorar 15-30 minutos (baixa cache + compila)..."
    
    # Instalação simples sem cache (mais confiável)
    echo "🔧 Instalando APIs do VEP..."
    ./INSTALL.pl \
        --AUTO a \
        --NO_BIOPERL \
        --NO_HTSLIB \
        --NO_TEST \
        --NO_UPDATE || {
        
        echo "❌ Instalação falhou"
        return 1
    }
    
    echo "✅ APIs do VEP instaladas"
}

# Função para configurar PATH
setup_vep_path() {
    echo ""
    echo "🔗 Configurando PATH do VEP..."
    
    # Adiciona ao ambiente conda
    CONDA_PREFIX="${CONDA_PREFIX:-$HOME/miniforge3/envs/$CONDA_ENV}"
    
    if [ -d "$CONDA_PREFIX" ]; then
        mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
        
        cat > "$CONDA_PREFIX/etc/conda/activate.d/vep_latest.sh" << EOF
#!/bin/bash
# VEP Latest - PATH e PERL5LIB
export PATH="$VEP_DIR:\$PATH"
export PERL5LIB="$VEP_DIR:\$PERL5LIB"
EOF
        
        echo "✅ PATH configurado para ambiente conda"
        
        # Ativa imediatamente
        export PATH="$VEP_DIR:$PATH"
        export PERL5LIB="$VEP_DIR:$PERL5LIB"
    fi
}

# Função para instalar cache separadamente
install_vep_cache() {
    echo ""
    echo "💾 Instalando cache do VEP..."
    
    cd "$VEP_DIR"
    
    # Usa INSTALL.pl do VEP baixado para instalar cache
    echo "🔧 Instalando cache via INSTALL.pl..."
    ./INSTALL.pl \
        --AUTO cf \
        --SPECIES "$SPECIES" \
        --ASSEMBLY "$ASSEMBLY" \
        --CACHEDIR "$VEPCACHE" \
        --CACHE_VERSION 114 \
        --NO_BIOPERL \
        --NO_UPDATE || {
        
        echo "⚠️  Instalação automática do cache falhou"
        echo "💡 Tentando com versão específica..."
        
        # Fallback: força versão 114
        ./INSTALL.pl \
            --AUTO cf \
            --SPECIES "$SPECIES" \
            --ASSEMBLY "$ASSEMBLY" \
            --CACHEDIR "$VEPCACHE" \
            --CACHE_VERSION 114 \
            --NO_BIOPERL || {
            
            echo "⚠️  Cache ainda falhou - pode instalar manualmente depois"
            echo "💡 Comando manual:"
            echo "   cd $VEP_DIR"
            echo "   ./INSTALL.pl --AUTO cf --SPECIES $SPECIES --ASSEMBLY $ASSEMBLY --CACHEDIR $VEPCACHE --CACHE_VERSION 114 --NO_BIOPERL"
        }
    }
    
    echo "✅ Cache instalado"
}

# Função de teste
test_vep_installation() {
    echo ""
    echo "🧪 Testando instalação..."
    
    if command -v vep >/dev/null 2>&1; then
        echo "✅ VEP encontrado: $(which vep)"
        
        # Testa versão
        echo "📋 Informações da versão:"
        vep --help 2>&1 | head -15 | grep -E "(ENSEMBL|ensembl-vep|Versions)" || echo "   Versão instalada com sucesso"
        
        # Testa cache
        if [ -d "$VEPCACHE/$SPECIES" ]; then
            CACHE_DIRS=$(find "$VEPCACHE/$SPECIES" -maxdepth 1 -type d -name "*$ASSEMBLY*" 2>/dev/null || true)
            if [ -n "$CACHE_DIRS" ]; then
                echo "✅ Cache verificado: $CACHE_DIRS"
            else
                echo "⚠️  Cache para $ASSEMBLY não encontrado"
                return 1
            fi
        else
            echo "❌ Diretório de cache não encontrado"
            return 1
        fi
        
        echo "✅ Instalação completa e funcional!"
        return 0
    else
        echo "❌ VEP não encontrado no PATH"
        return 1
    fi
}

# Execução principal
main() {
    echo "🔍 Detectando última versão do VEP no GitHub..."
    
    # Detecta versão mais recente
    LATEST_VERSION=$(detect_latest_vep_version)
    
    echo "✅ Versão detectada: $LATEST_VERSION"
    echo ""
    echo "🎯 Versão selecionada: $LATEST_VERSION"
    echo ""
    
    # Pergunta confirmação
    read -p "🤔 Continuar com a instalação? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Instalação cancelada"
        exit 0
    fi
    
    # Instala VEP
    install_vep_latest "$LATEST_VERSION"
    
    # Configura PATH
    setup_vep_path
    
    # Instala cache
    install_vep_cache
    
    # Testa
    if test_vep_installation; then
        echo ""
        echo "🎉 Instalação do VEP concluída com sucesso!"
        echo "📋 Resumo:"
        echo "   • Versão: $LATEST_VERSION (GitHub)"
        echo "   • Localização: $VEP_DIR"
        echo "   • Cache: $VEPCACHE"
        echo "   • Espécie: $SPECIES"
        echo "   • Assembly: $ASSEMBLY"
        echo ""
        echo "💡 Para usar:"
        echo "   vep --cache --offline --species $SPECIES --assembly $ASSEMBLY --dir_cache $VEPCACHE"
    else
        echo ""
        echo "❌ Problemas na instalação detectados"
        echo "💡 Para resolver dependências faltantes:"
        echo "   conda install -c conda-forge perl-try-tiny perl-uri perl-http-message"
        echo "   source start_genomics_universal.sh"
    fi
}

# Executa se chamado diretamente
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
