#!/bin/bash
# fix_vep_dependencies.sh - Corrige dependências Perl do VEP

echo "🔧 Correção Rápida - Dependências VEP"
echo "====================================="

# Verifica ambiente
if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "❌ Ambiente conda não ativo"
    echo "💡 Execute: source start_genomics_universal.sh"
    exit 1
fi

echo "✅ Ambiente ativo: $CONDA_DEFAULT_ENV"

echo ""
echo "📦 Instalando dependências Perl que estavam faltando..."

# Lista das dependências que causaram erro
MISSING_DEPS=(
    "perl-try-tiny"           # Try::Tiny
    "perl-uri"                # URI
    "perl-http-message"       # HTTP::Message  
    "perl-io-string"          # IO::String
    "perl-text-csv"           # Text::CSV
    "perl-file-temp"          # File::Temp
    "perl-getopt-long"        # Getopt::Long
)

echo "🔧 Instalando dependências críticas..."
for dep in "${MISSING_DEPS[@]}"; do
    echo "   📦 $dep"
    conda install -c conda-forge -c bioconda "$dep" -y -q || echo "   ⚠️  $dep falhou"
done

echo ""
echo "🧪 Testando VEP após correções..."

# Testa se VEP funciona agora
VEP_DIR="$HOME/ensembl-vep"
if [ -d "$VEP_DIR" ]; then
    export PATH="$VEP_DIR:$PATH"
    export PERL5LIB="$VEP_DIR:$PERL5LIB"
fi

if command -v vep >/dev/null 2>&1; then
    echo "✅ VEP encontrado: $(which vep)"
    
    # Testa help básico
    if vep --help >/dev/null 2>&1; then
        echo "✅ VEP funcionando!"
        
        # Mostra versão
        echo "📋 Versão:"
        vep 2>&1 | head -10 | grep -E "ensembl-vep|Versions" || echo "   VEP instalado com sucesso"
        
    else
        echo "❌ VEP ainda com problemas"
        echo "🔍 Testando dependências específicas..."
        
        # Testa dependências uma por uma
        PERL_MODULES=("Try::Tiny" "List::MoreUtils" "URI" "HTTP::Message")
        for module in "${PERL_MODULES[@]}"; do
            if perl -M"$module" -e 1 2>/dev/null; then
                echo "   ✅ $module"
            else
                echo "   ❌ $module - FALTANDO"
            fi
        done
    fi
else
    echo "❌ VEP não encontrado"
fi

echo ""
echo "💡 Se ainda houver problemas:"
echo "   1. Reinstale VEP via conda: conda install -c bioconda ensembl-vep --force-reinstall"
echo "   2. Use instalador original: source vep_install.sh"
echo "   3. Instale dependência crítica: conda install -c conda-forge perl-try-tiny"

echo ""
echo "🎯 Para continuar com pipeline:"
echo "   source start_genomics_universal.sh"
echo "   ./genomes_analyzer.py --config config_emergency_fix.yaml"

