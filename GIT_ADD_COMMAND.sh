#!/bin/bash
# Comando Git Add para Neural Module - Versão Organizada
# Gerado automaticamente em: Outubro 2025

echo "📦 Preparando arquivos do Neural Module para commit..."
echo ""

# =====================================================
# 📚 DOCUMENTAÇÃO PRINCIPAL (NOVA ESTRUTURA)
# =====================================================
echo "1. Documentação Principal..."
git add \
  NEURAL_MODULE.md \
  INSTALL_NEURAL.md \
  USAGE_NEURAL.md \
  RESULTS_NEURAL.md \
  NEURAL_CHANGELOG.md

# =====================================================
# 📖 DOCUMENTAÇÃO EXISTENTE
# =====================================================
echo "2. Documentação Técnica..."
git add \
  NEURAL_MODULE_README.md \
  NEURAL_QUICKSTART.md \
  NEURAL_MODULE_INDEX.md \
  LEIA-ME_NEURAL.md \
  OUTPUTS_DISPONIVEIS.md \
  TAMANHOS_SUPORTADOS.md \
  CORRECOES_APLICADAS.md \
  VISUALIZACAO_CORRIGIDA.md \
  VISUALIZACOES_AVANCADAS.md

# =====================================================
# 🐍 MÓDULOS PYTHON
# =====================================================
echo "3. Módulos Python..."
git add \
  neural_module.py \
  neural_example.py \
  neural_integration.py \
  neural_visualizations_advanced.py

# =====================================================
# 🧬 EXEMPLO BIOLÓGICO
# =====================================================
echo "4. Exemplo de Sequência (Gene HBB)..."
git add example_sequence.fasta

# =====================================================
# ⚙️ CONFIGURAÇÃO
# =====================================================
echo "5. Configuração..."
git add neural_config.yaml

# =====================================================
# 🔧 SCRIPTS DE INSTALAÇÃO E TESTE
# =====================================================
echo "6. Scripts de Instalação..."
git add \
  install_alphagenome.sh \
  check_neural_requirements.sh

echo "7. Scripts de Teste e Demo..."
git add \
  test_neural_module.sh \
  demo_neural_module.sh \
  show_neural_summary.sh

# =====================================================
# 🔍 SCRIPTS DE DIAGNÓSTICO
# =====================================================
echo "8. Scripts de Diagnóstico..."
git add \
  check_alphagenome_outputs.py \
  check_dna_client_methods.py \
  check_visualization_api.py \
  check_output_structure.py

# =====================================================
# 📄 ESTE SCRIPT
# =====================================================
echo "9. Script de Git Add..."
git add GIT_ADD_COMMAND.sh

echo ""
echo "✅ Todos os arquivos adicionados ao staging!"
echo ""
echo "📊 Resumo:"
echo "   • Documentação Principal: 5 arquivos"
echo "   • Documentação Técnica: 9 arquivos"
echo "   • Módulos Python: 4 arquivos"
echo "   • Exemplo Biológico: 1 arquivo"
echo "   • Configuração: 1 arquivo"
echo "   • Scripts: 9 arquivos"
echo "   • Total: 29 arquivos"
echo ""
echo "📝 Próximo passo:"
echo "   git status  # Verificar arquivos staged"
echo ""
echo "💡 Sugestão de commit:"
echo '   git commit -m "feat: Reorganize neural_module documentation and enhance defaults'
echo ''
echo '   - NEW: Structured documentation (NEURAL_MODULE.md + 3 guides)'
echo '   - NEW: Real biological example (HBB gene, sickle cell region)'
echo '   - CHANGED: Advanced visualizations now default'
echo '   - CHANGED: Ontology information displayed automatically'
echo '   - IMPROVED: Simplified CLI (removed --advanced-viz flag)'
echo '   - IMPROVED: Better user experience with automatic features'
echo ''
echo '   Complete documentation restructure with installation, usage,'
echo '   and results interpretation guides. Example sequence now contains'
echo '   real HBB gene region for sickle cell anemia variant analysis.'
echo '   Advanced visualizations and ontology information enabled by default.'
echo ''
echo '   Total: 29 files | ~4,500 lines of code and documentation"'
echo ""

