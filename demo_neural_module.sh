#!/usr/bin/env bash
# Demonstração rápida do Neural Module

cat << 'EOF'
╔═══════════════════════════════════════════════════════════════════╗
║         Neural Module - Demonstração de Uso                      ║
║         Análise de DNA com AlphaGenome (Google DeepMind)         ║
╚═══════════════════════════════════════════════════════════════════╝

Este script demonstra os principais recursos do neural_module.py:

1️⃣  Análise Básica de Sequências
   Analisa sequências de DNA e prediz características funcionais

2️⃣  Análise com Múltiplos Outputs
   RNA-seq, ATAC-seq, marcadores de histonas, etc.

3️⃣  Análise de Variantes
   Prediz o efeito de variantes SNP (A>C, G>T, etc.)

4️⃣  Visualizações Avançadas
   Gera gráficos em alta resolução (PNG, PDF, SVG)

───────────────────────────────────────────────────────────────────

📋 PRÉ-REQUISITOS:

1. Instalar AlphaGenome:
   $ bash install_alphagenome.sh

2. Obter API key (gratuita para uso não comercial):
   👉 https://www.alphagenomedocs.com/

3. Preparar arquivo FASTA com suas sequências

───────────────────────────────────────────────────────────────────

🚀 EXEMPLOS DE USO:

# Exemplo 1: Análise básica
$ python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/

# Exemplo 2: Análise específica de RNA-seq e ATAC-seq
$ python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \\
    --outputs RNA_SEQ ATAC

# Exemplo 3: Análise de variante A→C na posição 1000
$ python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \\
    --variant 1000 A C

# Exemplo 4: Alta resolução com múltiplos formatos
$ python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \\
    --formats png pdf svg --dpi 600

# Exemplo 5: Análise completa de região genômica
$ python neural_module.py -i chr1_region.fasta -k YOUR_API_KEY -o chr1/ \\
    --chromosome chr1 --start 1000000 \\
    --outputs RNA_SEQ CAGE ATAC H3K27AC H3K4ME3 H3K27ME3 CTCF \\
    --formats png pdf --dpi 600

───────────────────────────────────────────────────────────────────

📊 OUTPUTS DISPONÍVEIS:

Expressão Gênica:
  • RNA_SEQ  - Predição via RNA-seq
  • CAGE     - Cap Analysis of Gene Expression

Cromatina:
  • ATAC     - Acessibilidade de cromatina

Marcadores de Histonas:
  • H3K27AC  - Elementos regulatórios ativos (enhancers)
  • H3K4ME3  - Promotores ativos
  • H3K27ME3 - Repressão gênica (Polycomb)
  • H3K36ME3 - Corpos gênicos ativos
  • H3K9ME3  - Heterocromatina

Fatores de Transcrição:
  • CTCF     - Fator de ligação de insulador

───────────────────────────────────────────────────────────────────

📁 ESTRUTURA DE SAÍDA:

results/
├── sequence_id_1_RNA_SEQ.png       # Gráfico de RNA-seq
├── sequence_id_1_ATAC.png          # Gráfico de ATAC-seq
├── sequence_id_1_H3K27AC.png       # Gráfico de H3K27AC
├── sequence_id_2_RNA_SEQ.png
├── ...
└── analysis_report.json            # Relatório da análise

───────────────────────────────────────────────────────────────────

💡 CASOS DE USO:

1. Pesquisa de Variantes Regulatórias
   Identificar impacto de SNPs em regiões não codificantes

2. Predição de Elementos Funcionais
   Identificar promotores, enhancers, insuladores

3. Estudos de Expressão Diferencial
   Comparar predições entre variantes

4. Análise de Cromatina
   Estudar acessibilidade e modificações epigenéticas

5. Genômica Funcional
   Caracterizar sequências de DNA desconhecidas

───────────────────────────────────────────────────────────────────

⚙️  OPÇÕES AVANÇADAS:

--chromosome CHR      Cromossomo de referência (padrão: chr1)
--start POS           Posição inicial (padrão: 1000000)
--variant POS REF ALT Analisar variante específica
--formats png pdf svg Formatos de saída dos gráficos
--dpi 600             Resolução dos gráficos
--no-plots            Apenas análise, sem gráficos

───────────────────────────────────────────────────────────────────

📚 RECURSOS:

Documentação:  https://www.alphagenomedocs.com/
GitHub:        https://github.com/google-deepmind/alphagenome
Paper:         Avsec et al. 2025 - "AlphaGenome: advancing 
               regulatory variant effect prediction"

───────────────────────────────────────────────────────────────────

🤝 SUPORTE:

README:        NEURAL_MODULE_README.md
Testes:        bash test_neural_module.sh YOUR_API_KEY
Exemplo FASTA: example_sequence.fasta

───────────────────────────────────────────────────────────────────

⚠️  NOTAS IMPORTANTES:

• Uso gratuito apenas para pesquisa não comercial
• Requer conexão com internet (API online)
• Taxa de queries limitada (adequado para ~1000s de predições)
• Sequências: 100 bp a 1 Mbp (1 milhão de pares de bases)
• Nunca compartilhe sua API key publicamente

───────────────────────────────────────────────────────────────────

🎯 TESTE RÁPIDO:

Para testar rapidamente com a sequência de exemplo:

$ python neural_module.py \\
    -i example_sequence.fasta \\
    -k YOUR_API_KEY \\
    -o test_results/ \\
    --outputs RNA_SEQ ATAC

╔═══════════════════════════════════════════════════════════════════╗
║  Pronto para começar! Substitua YOUR_API_KEY e execute! 🚀       ║
╚═══════════════════════════════════════════════════════════════════╝

EOF

