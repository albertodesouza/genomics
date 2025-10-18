#!/bin/bash
# vep_update_to_latest.sh - Atualiza VEP para última versão (método híbrido)

echo "🔄 Atualizador VEP - Método Híbrido (Conda + GitHub)"
echo "=================================================="

SPECIES="homo_sapiens"
ASSEMBLY="GRCh38"
VEPCACHE="/dados/vep_cache"

echo "📋 Estratégia:"
echo "   1. Usar VEP do conda (mais estável)"
echo "   2. Atualizar cache para versão mais recente"
echo "   3. Configurar para usar cache atualizado"

# Verifica ambiente
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "❌ Ambiente conda não ativo"
    echo "💡 Execute: source start_genomics_universal.sh"
    exit 1
fi

echo ""
echo "🔧 1. Atualizando VEP via conda..."

# Atualiza VEP e dependências
conda update -c bioconda ensembl-vep -y || echo "⚠️  Update falhou, continuando..."

# Instala dependências que faltam
conda install -c bioconda -c conda-forge -y \
    perl-list-moreutils \
    perl-set-intervaltree \
    perl-bio-db-hts || echo "⚠️  Algumas dependências falharam"

echo ""
echo "🧪 2. Testando VEP atual..."

if command -v vep >/dev/null 2>&1; then
    echo "✅ VEP encontrado: $(which vep)"
    
    # Testa funcionalidade básica
    if vep --help >/dev/null 2>&1; then
        echo "✅ VEP funcional"
        VEP_INFO=$(vep 2>&1 | head -15 | grep -E "(ensembl-vep|Versions)" | head -3 || echo "VEP conda version")
        echo "📋 Versão atual:"
        echo "$VEP_INFO"
    else
        echo "❌ VEP com problemas"
        echo "💡 Instalando dependência crítica..."
        conda install -c conda-forge perl-list-moreutils -y
    fi
else
    echo "❌ VEP não encontrado"
    echo "🔧 Instalando VEP via conda..."
    conda install -c bioconda ensembl-vep -y
fi

echo ""
echo "💾 3. Configurando cache atualizado..."

mkdir -p "$VEPCACHE"

# Tenta baixar cache mais recente
if command -v vep >/dev/null 2>&1; then
    echo "📥 Baixando cache para $SPECIES $ASSEMBLY..."
    
    # Método preferido: vep_install
    if command -v vep_install >/dev/null 2>&1; then
        echo "🎯 Usando vep_install..."
        vep_install -a cf -s "$SPECIES" -y "$ASSEMBLY" -c "$VEPCACHE" --NO_BIOPERL --CONVERT || {
            echo "⚠️  vep_install falhou, cache pode estar desatualizado"
        }
    else
        echo "⚠️  vep_install não disponível"
        echo "💡 Cache manual necessário - use o instalador original se precisar"
    fi
fi

echo ""
echo "🧪 4. Verificação final..."

# Teste completo
if command -v vep >/dev/null 2>&1 && vep --help >/dev/null 2>&1; then
    echo "✅ VEP funcionando"
    
    # Verifica cache
    if [ -d "$VEPCACHE/$SPECIES" ]; then
        CACHE_DIRS=$(find "$VEPCACHE/$SPECIES" -maxdepth 1 -type d -name "*$ASSEMBLY*" 2>/dev/null)
        if [ -n "$CACHE_DIRS" ]; then
            echo "✅ Cache encontrado:"
            echo "$CACHE_DIRS"
            
            # Mostra versão do cache
            for cache_dir in $CACHE_DIRS; do
                if [ -f "$cache_dir/info.txt" ]; then
                    echo "📋 Info do cache:"
                    head -5 "$cache_dir/info.txt" 2>/dev/null || echo "   Cache válido"
                fi
            done
        else
            echo "⚠️  Cache para $ASSEMBLY não encontrado"
        fi
    else
        echo "⚠️  Diretório de cache não existe"
    fi
    
    echo ""
    echo "🎯 VEP Atualizado!"
    echo "================"
    echo "✅ Status: Funcional"
    echo "📁 Cache: $VEPCACHE"
    echo "🧬 Espécie: $SPECIES"
    echo "🧬 Assembly: $ASSEMBLY"
    echo ""
    echo "💡 Para usar no pipeline:"
    echo "   vep --cache --offline --species $SPECIES --assembly $ASSEMBLY --dir_cache $VEPCACHE"
    
else
    echo "❌ VEP ainda com problemas"
    echo ""
    echo "🛠️  Soluções:"
    echo "   1. Reinstalar dependências: conda install -c conda-forge perl-list-moreutils"
    echo "   2. Usar instalador original: source vep_install.sh"
    echo "   3. Verificar ambiente: source start_genomics_universal.sh"
fi
