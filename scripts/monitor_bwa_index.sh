#!/bin/bash
# monitor_bwa_index.sh - Monitora progresso do BWA index em tempo real

echo "📊 Monitor de Progresso BWA Index"
echo "================================"

# Detecta qual processo BWA está rodando
BWA_PID=$(pgrep -f "bwa.*index" | head -1)
if [ -z "$BWA_PID" ]; then
    echo "❌ Nenhum processo BWA index encontrado"
    exit 1
fi

BWA_CMD=$(ps -p "$BWA_PID" -o cmd --no-headers)
echo "🔍 Processo detectado (PID: $BWA_PID):"
echo "   $BWA_CMD"
echo ""

# Diretório de referência
REF_DIR="/mnt/barra-dados/GENOMICS_DATA/monster_256gb/refs"
if [ ! -d "$REF_DIR" ]; then
    REF_DIR="refs"
fi

echo "📁 Monitorando arquivos em: $REF_DIR"
echo "⏱️  Atualizações a cada 30 segundos (Ctrl+C para parar)"
echo ""

START_TIME=$(date +%s)

while kill -0 "$BWA_PID" 2>/dev/null; do
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))
    ELAPSED_MIN=$((ELAPSED / 60))
    ELAPSED_SEC=$((ELAPSED % 60))
    
    # Verifica arquivos sendo criados
    FILES_INFO=()
    TOTAL_SIZE=0
    
    # BWA clássico
    for ext in amb ann bwt pac sa; do
        FILE="$REF_DIR/reference.fa.$ext"
        if [ -f "$FILE" ]; then
            SIZE=$(stat -c%s "$FILE" 2>/dev/null || echo 0)
            TOTAL_SIZE=$((TOTAL_SIZE + SIZE))
            SIZE_HUMAN=$(numfmt --to=iec-i --suffix=B "$SIZE" 2>/dev/null || echo "${SIZE}B")
            FILES_INFO+=("$ext($SIZE_HUMAN)")
        fi
    done
    
    # BWA-MEM2
    for ext in 0123 bwt.2bit.64; do
        FILE="$REF_DIR/reference.fa.$ext"
        if [ -f "$FILE" ]; then
            SIZE=$(stat -c%s "$FILE" 2>/dev/null || echo 0)
            TOTAL_SIZE=$((TOTAL_SIZE + SIZE))
            SIZE_HUMAN=$(numfmt --to=iec-i --suffix=B "$SIZE" 2>/dev/null || echo "${SIZE}B")
            EXT_DISPLAY=$(echo "$ext" | sed 's/\.2bit\.64/2bit/')
            FILES_INFO+=("$EXT_DISPLAY($SIZE_HUMAN)")
        fi
    done
    
    # Mostra progresso
    TOTAL_HUMAN=$(numfmt --to=iec-i --suffix=B "$TOTAL_SIZE" 2>/dev/null || echo "${TOTAL_SIZE}B")
    
    if [ ${#FILES_INFO[@]} -gt 0 ]; then
        STATUS="total $TOTAL_HUMAN • $(IFS=', '; echo "${FILES_INFO[*]:0:3}")"
        if [ ${#FILES_INFO[@]} -gt 3 ]; then
            STATUS="$STATUS +$((${#FILES_INFO[@]} - 3)) mais"
        fi
    else
        STATUS="iniciando..."
    fi
    
    printf "\r[$(date '+%H:%M:%S')] BWA index … %dm%02ds • %s" \
           "$ELAPSED_MIN" "$ELAPSED_SEC" "$STATUS"
    
    sleep 30
done

echo ""
echo ""
echo "✅ Processo BWA index concluído!"

# Relatório final
if [ -d "$REF_DIR" ]; then
    echo "📊 Arquivos finais criados:"
    ls -lh "$REF_DIR"/reference.fa.* 2>/dev/null | while read line; do
        echo "   $line"
    done
    
    FINAL_SIZE=$(du -sh "$REF_DIR"/reference.fa.* 2>/dev/null | awk '{sum+=$1} END {print sum}' || echo "0")
    echo ""
    echo "💾 Tamanho total dos índices: $(du -ch "$REF_DIR"/reference.fa.* 2>/dev/null | tail -1 | cut -f1)"
fi
