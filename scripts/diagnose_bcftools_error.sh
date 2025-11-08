#!/bin/bash
# diagnose_bcftools_error.sh - Diagnostica erros específicos do bcftools

echo "🔍 Diagnóstico de Erros BCFtools"
echo "================================"

# Encontra logs de erro mais recentes
echo "📝 Procurando logs de erro..."
ERROR_LOGS=$(find logs/ -name "bcftools_*.log" -newer logs/bcftools_NA12878_part_01.log 2>/dev/null | head -5)

if [ -z "$ERROR_LOGS" ]; then
    echo "❌ Nenhum log de erro encontrado"
    echo "💡 Verifique: ls -la logs/bcftools_*"
    exit 1
fi

echo "📋 Logs encontrados:"
for log in $ERROR_LOGS; do
    echo "   📄 $log"
done

echo ""
echo "🔍 Analisando erros..."

for log in $ERROR_LOGS; do
    echo ""
    echo "📄 Arquivo: $log"
    echo "----------------------------------------"
    
    # Verifica se há erros específicos
    if grep -q "error\|Error\|ERROR\|failed\|Failed\|FAILED" "$log"; then
        echo "❌ ERROS ENCONTRADOS:"
        grep -i "error\|failed" "$log" | head -10
    fi
    
    # Verifica cromossomos problemáticos
    if grep -q "EBV\|CMV\|HIV\|HPV\|HCV\|HTLV\|KSHV" "$log"; then
        echo "🦠 CROMOSSOMOS VIRAIS DETECTADOS:"
        grep "EBV\|CMV\|HIV\|HPV\|HCV\|HTLV\|KSHV" "$log" | head -5
    fi
    
    # Mostra últimas linhas
    echo "📖 Últimas linhas:"
    tail -10 "$log"
    echo "----------------------------------------"
done

echo ""
echo "🔍 Verificando shard problemático..."

# Identifica qual shard falhou
FAILED_SHARD=$(echo "$ERROR_LOGS" | head -1 | sed 's/.*part_\([0-9]*\)\.log/\1/')
if [ -n "$FAILED_SHARD" ]; then
    BED_FILE="vcf/shards/NA12878/part_${FAILED_SHARD}.bed"
    if [ -f "$BED_FILE" ]; then
        echo "📋 Cromossomos no shard $FAILED_SHARD:"
        echo "----------------------------------------"
        cut -f1 "$BED_FILE" | sort -u | head -20
        echo "----------------------------------------"
        
        echo ""
        echo "🦠 Cromossomos virais no shard:"
        cut -f1 "$BED_FILE" | sort -u | grep -E "EBV|CMV|HIV|HPV|HCV|HTLV|KSHV" || echo "   (nenhum detectado)"
    fi
fi

echo ""
echo "💡 Soluções sugeridas:"
echo "   1. Adicionar filter_problematic_contigs: true no YAML"
echo "   2. Usar config_human_30x_filtered.yaml"
echo "   3. Verificar logs específicos acima para detalhes"
