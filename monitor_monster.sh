#!/bin/bash
# monitor_monster.sh - Script para monitorar pipeline em background

PID_FILE="pipeline_monster.pid"

if [ ! -f "$PID_FILE" ]; then
    echo "❌ Pipeline não está rodando (PID file não encontrado)"
    exit 1
fi

PID=$(cat "$PID_FILE")

if ! kill -0 "$PID" 2>/dev/null; then
    echo "❌ Pipeline não está mais rodando (PID $PID morto)"
    rm -f "$PID_FILE"
    exit 1
fi

# Encontra o log mais recente
LOG_FILE=$(ls -t pipeline_monster_*.log 2>/dev/null | head -1)

if [ -z "$LOG_FILE" ]; then
    echo "❌ Arquivo de log não encontrado"
    exit 1
fi

echo "🔍 Monitorando pipeline monster (PID: $PID)"
echo "📝 Log: $LOG_FILE"
echo "📊 Uso de recursos:"

# Mostra uso de CPU e memória
ps -p "$PID" -o pid,ppid,pcpu,pmem,etime,cmd --no-headers 2>/dev/null || echo "Processo não encontrado"

echo ""
echo "📖 Últimas 20 linhas do log:"
echo "----------------------------------------"
tail -20 "$LOG_FILE"
echo "----------------------------------------"
echo ""
echo "💡 Para monitoramento contínuo:"
echo "   tail -f $LOG_FILE"
echo "   watch 'tail -20 $LOG_FILE'"
echo "   htop -p $PID"
