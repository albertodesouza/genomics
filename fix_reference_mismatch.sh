#!/bin/bash
# fix_reference_mismatch.sh - Corrige incompatibilidades de referência

set -euo pipefail

echo "🔧 Diagnóstico e Correção de Reference Mismatch"
echo "=============================================="

# Verifica arquivos atuais
echo "📁 Verificando arquivos atuais..."
REF_FASTA="/mnt/barra-dados/GENOMICS_DATA/low_memory/refs/reference.fa"
BAM_FILE="/mnt/barra-dados/GENOMICS_DATA/low_memory/bam/NA12878.mkdup.bam"

if [ ! -f "$REF_FASTA" ]; then
    echo "❌ FASTA de referência não encontrado: $REF_FASTA"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "❌ BAM não encontrado: $BAM_FILE"
    exit 1
fi

echo "✅ Arquivos encontrados"

# Verifica headers do BAM para ver qual referência foi usada
echo ""
echo "🔍 Verificando header do BAM..."
samtools view -H "$BAM_FILE" | grep "^@SQ" | head -5
echo ""

# Verifica se há incompatibilidade na posição específica do erro
echo "🔍 Verificando posição específica do erro (chr22:45756861)..."
echo "Referência FASTA:"
samtools faidx "$REF_FASTA" chr22:45756861-45756861
echo ""

echo "Reads no BAM nesta posição:"
samtools view "$BAM_FILE" chr22:45756861-45756861 | head -3 || echo "Nenhum read encontrado"
echo ""

# Opções de solução
echo "💡 Soluções disponíveis:"
echo ""
echo "1️⃣  SOLUÇÃO RÁPIDA - Usar configuração com referência atualizada:"
echo "   ./genomes_analyzer.py --config config_human_30x_latest_ref.yaml"
echo "   ⚠️  Isso irá:"
echo "   • Baixar referência GENCODE release 46 (mais recente)"
echo "   • Recriar índices BWA"
echo "   • Realinhar todas as amostras (pode demorar)"
echo ""
echo "2️⃣  SOLUÇÃO CONSERVADORA - Usar apenas cromossomos principais:"
echo "   ./genomes_analyzer.py --config config_emergency_fix.yaml"
echo "   ✅ Isso irá:"
echo "   • Usar apenas chr1-22,X,Y,MT"
echo "   • Evitar regiões problemáticas"
echo "   • Manter alinhamentos existentes"
echo ""
echo "3️⃣  SOLUÇÃO MANUAL - Corrigir posição específica:"
echo "   # Editar FASTA para corrigir posição chr22:45756861"
echo "   # (não recomendado - pode causar outros problemas)"
echo ""

# Recomendação
echo "🎯 RECOMENDAÇÃO:"
echo "   Use a Solução 1 (config_human_30x_latest_ref.yaml) para:"
echo "   • Referência atualizada e compatível"
echo "   • Zero problemas de mismatch"
echo "   • Resultados mais confiáveis"
echo ""
echo "   ⏱️  Tempo extra: ~2-4h para realinhamento"
echo "   💾 Espaço extra: ~50GB para nova referência"
echo ""
echo "   Ou use Solução 2 (config_emergency_fix.yaml) para:"
echo "   • Solução rápida com dados existentes"
echo "   • Pode perder algumas variantes em regiões problemáticas"
