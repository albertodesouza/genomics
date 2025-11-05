#!/usr/bin/env bash
# Script de diagnóstico para verificar requisitos do neural_module.py

set -e

ERRORS=0
WARNINGS=0

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  Verificação de Requisitos - Neural Module                       ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

# Função para verificar comando
check_command() {
    local cmd=$1
    local name=$2
    local required=$3
    
    if command -v "$cmd" &> /dev/null; then
        echo "✅ $name: instalado"
        if [ "$cmd" = "python" ] || [ "$cmd" = "python3" ]; then
            local version=$($cmd --version 2>&1)
            echo "   Versão: $version"
        fi
        return 0
    else
        if [ "$required" = "yes" ]; then
            echo "❌ $name: NÃO INSTALADO (obrigatório)"
            ((ERRORS++))
        else
            echo "⚠️  $name: não instalado (opcional)"
            ((WARNINGS++))
        fi
        return 1
    fi
}

# Função para verificar módulo Python
check_python_module() {
    local module=$1
    local name=$2
    local required=$3
    
    if python3 -c "import $module" 2>/dev/null; then
        echo "✅ $name: instalado"
        # Tentar pegar versão
        local version=$(python3 -c "import $module; print(getattr($module, '__version__', 'N/A'))" 2>/dev/null || echo "N/A")
        if [ "$version" != "N/A" ]; then
            echo "   Versão: $version"
        fi
        return 0
    else
        if [ "$required" = "yes" ]; then
            echo "❌ $name: NÃO INSTALADO (obrigatório)"
            echo "   Instalar: pip install $module"
            ((ERRORS++))
        else
            echo "⚠️  $name: não instalado (opcional)"
            echo "   Instalar: pip install $module"
            ((WARNINGS++))
        fi
        return 1
    fi
}

# Verificar arquivos principais
echo "📁 Verificando arquivos do projeto..."
echo ""

files=("neural_module.py" "neural_example.py" "neural_integration.py" "example_sequence.fasta")
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ $file: encontrado"
    else
        echo "❌ $file: NÃO ENCONTRADO"
        ((ERRORS++))
    fi
done

echo ""
echo "🐍 Verificando Python e módulos básicos..."
echo ""

# Python
check_command python3 "Python 3" yes || check_command python "Python" yes

# Módulos Python básicos
check_python_module rich "rich" yes
check_python_module matplotlib "matplotlib" yes
check_python_module numpy "numpy" yes
check_python_module pandas "pandas" no

echo ""
echo "🧬 Verificando AlphaGenome..."
echo ""

# AlphaGenome
if check_python_module alphagenome "alphagenome" yes; then
    echo "   ℹ️  AlphaGenome está instalado!"
    echo ""
    echo "   Componentes do AlphaGenome:"
    
    # Verificar submodulos
    if python3 -c "from alphagenome.models import dna_client" 2>/dev/null; then
        echo "   ✅ alphagenome.models.dna_client"
    else
        echo "   ❌ alphagenome.models.dna_client"
        ((ERRORS++))
    fi
    
    if python3 -c "from alphagenome.data import genome" 2>/dev/null; then
        echo "   ✅ alphagenome.data.genome"
    else
        echo "   ❌ alphagenome.data.genome"
        ((ERRORS++))
    fi
    
    if python3 -c "from alphagenome.visualization import plot_components" 2>/dev/null; then
        echo "   ✅ alphagenome.visualization.plot_components"
    else
        echo "   ⚠️  alphagenome.visualization.plot_components (opcional)"
        ((WARNINGS++))
    fi
else
    echo "   ❌ AlphaGenome não está instalado!"
    echo ""
    echo "   Para instalar, execute:"
    echo "   $ bash install_alphagenome.sh"
    echo "   ou"
    echo "   $ git clone https://github.com/google-deepmind/alphagenome.git"
    echo "   $ pip install ./alphagenome"
fi

echo ""
echo "🔧 Verificando ferramentas opcionais..."
echo ""

# Ferramentas opcionais (para neural_integration.py)
check_command bcftools "bcftools" no
check_command samtools "samtools" no
check_command bedtools "bedtools" no

echo ""
echo "🌐 Verificando conectividade..."
echo ""

# Testar conectividade com Google (para API)
if ping -c 1 google.com &> /dev/null; then
    echo "✅ Conectividade com internet: OK"
else
    echo "⚠️  Conectividade com internet: FALHOU"
    echo "   A API do AlphaGenome requer conexão com internet"
    ((WARNINGS++))
fi

echo ""
echo "📊 Verificando ambiente conda..."
echo ""

if command -v conda &> /dev/null; then
    echo "✅ Conda: instalado"
    echo "   Versão: $(conda --version)"
    
    # Verificar ambiente genomics
    if conda env list | grep -q "^genomics "; then
        echo "✅ Ambiente 'genomics': encontrado"
    else
        echo "ℹ️  Ambiente 'genomics': não encontrado (opcional)"
    fi
else
    echo "ℹ️  Conda: não instalado (opcional para este módulo)"
fi

echo ""
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  RESUMO                                                           ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

if [ $ERRORS -eq 0 ]; then
    echo "✅ TODOS OS REQUISITOS OBRIGATÓRIOS FORAM ATENDIDOS!"
    echo ""
    
    if [ $WARNINGS -gt 0 ]; then
        echo "⚠️  $WARNINGS aviso(s) encontrado(s) (componentes opcionais)"
        echo ""
    fi
    
    echo "🚀 Você está pronto para usar o neural_module.py!"
    echo ""
    echo "Próximos passos:"
    echo "1. Obtenha sua API key em: https://www.alphagenomedocs.com/"
    echo "2. Execute: python neural_module.py -i example_sequence.fasta -k YOUR_API_KEY -o results/"
    echo ""
    exit 0
else
    echo "❌ $ERRORS ERRO(S) CRÍTICO(S) ENCONTRADO(S)"
    
    if [ $WARNINGS -gt 0 ]; then
        echo "⚠️  $WARNINGS aviso(s) encontrado(s)"
    fi
    
    echo ""
    echo "Por favor, corrija os erros acima antes de usar o neural_module.py"
    echo ""
    exit 1
fi

