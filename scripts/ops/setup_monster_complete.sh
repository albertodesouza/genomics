#!/bin/bash
# setup_monster_complete.sh - Instalação completa para máquina monster

set -euo pipefail

echo "🚀 Configuração Completa - Máquina Monster (128 cores + 256GB RAM)"
echo "=================================================================="

# Detecta usuário e cria diretórios
USER_HOME="$HOME"
echo "👤 Usuário: $(whoami)"
echo "🏠 Home: $USER_HOME"

# 1. Instala conda/mamba se não existir
echo ""
echo "📦 1. Verificando Conda/Mamba..."
if ! command -v conda >/dev/null 2>&1; then
    echo "🔧 Instalando conda/mamba..."
    chmod +x install_conda_universal.sh
    ./install_conda_universal.sh
    
    # Recarrega shell
    source ~/.bashrc || true
    export PATH="$USER_HOME/miniforge3/bin:$PATH"
else
    echo "✅ Conda já instalado: $(which conda)"
fi

# 2. Instala ambiente genomics
echo ""
echo "🧬 2. Configurando Ambiente Genomics..."
if ! conda env list | grep -q "^genomics "; then
    echo "🔧 Criando ambiente genomics..."
    chmod +x install_genomics_env.sh
    ./install_genomics_env.sh
else
    echo "✅ Ambiente genomics já existe"
fi

# 3. Configura VEP (versão mais recente)
echo ""
echo "🏷️  3. Configurando VEP (última versão do GitHub)..."
source start_genomics_universal.sh
if ! command -v vep >/dev/null 2>&1; then
    echo "🔧 Instalando VEP última versão..."
    chmod +x vep_install_smart.sh
    ./vep_install_smart.sh
else
    VEP_VERSION=$(vep --help 2>&1 | grep -E "ensembl-vep.*[0-9]" | head -1 || echo "versão desconhecida")
    echo "✅ VEP já instalado: $VEP_VERSION"
    echo "💡 Para atualizar para última versão: ./vep_install_smart.sh"
fi

# 4. Otimiza sistema para 256GB
echo ""
echo "⚙️  4. Otimizando Sistema para 256GB RAM..."
chmod +x setup_monster_256gb.sh
./setup_monster_256gb.sh

# 5. Cria diretórios de dados
echo ""
echo "📁 5. Criando Estrutura de Diretórios..."
sudo mkdir -p /dados/GENOMICS_DATA/monster_256gb
sudo mkdir -p /dados/vep_cache
sudo chown -R $(whoami):$(whoami) /dados/ || echo "⚠️  Aviso: não foi possível mudar owner de /dados"

# 6. Testa configuração
echo ""
echo "🧪 6. Testando Configuração..."
source start_genomics_universal.sh

echo "✅ Ferramentas verificadas:"
tools=("bwa-mem2" "samtools" "bcftools" "vep" "seqtk" "tabix")
all_ok=true
for tool in "${tools[@]}"; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "   ✅ $tool"
    else
        echo "   ❌ $tool - FALTANDO!"
        all_ok=false
    fi
done

# 7. Verifica recursos do sistema
echo ""
echo "💻 7. Recursos do Sistema:"
echo "   🧠 RAM: $(free -h | awk 'NR==2{print $2}')"
echo "   ⚡ CPU: $(nproc) cores"
echo "   💾 /dev/shm: $(df -h /dev/shm | awk 'NR==2{print $2}')"
echo "   📁 /dados: $(df -h /dados 2>/dev/null | awk 'NR==2{print $4}' || echo 'não montado')"

# 8. Instruções finais
echo ""
echo "🎯 Instalação Concluída!"
echo "======================="

if [ "$all_ok" = true ]; then
    echo "✅ Todos os componentes instalados com sucesso!"
    echo ""
    echo "🚀 Para executar pipeline:"
    echo "   ./run_monster_background.sh"
    echo ""
    echo "📊 Para monitorar:"
    echo "   ./monitor_monster.sh"
    echo ""
    echo "⚙️  Configuração otimizada para:"
    echo "   • 96 threads alinhamento"
    echo "   • 64 shards BCFtools paralelos"
    echo "   • 96 forks VEP"
    echo "   • 256GB RAM + RAM disk"
    echo "   • Filtro automático de cromossomos problemáticos"
else
    echo "❌ Alguns componentes faltando - verifique instalação"
    echo "💡 Tente executar novamente ou instale manualmente"
fi

echo ""
echo "📋 Arquivos de configuração criados:"
echo "   • config_human_30x_monster.yaml"
echo "   • run_monster_background.sh"
echo "   • monitor_monster.sh"
echo "   • start_genomics_universal.sh"
