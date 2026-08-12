#!/bin/bash
# run_in_background.sh - Script para executar pipeline em background

#set -euo pipefail

# Configurações
# Verificação se parâmetro foi fornecido
if [ -z "$1" ]; then
    echo "❌ Uso: $0 <arquivo_de_configuracao>"
    echo "   Exemplo: $0 config_human_30x_monster.yaml"
    exit 1
fi
CONFIG="$1"
LOG_FILE="pipeline_$(date +%Y%m%d_%H%M%S).log"
PID_FILE="pipeline.pid"

echo "🚀 Preparando execução em background"
echo "📋 Config: $CONFIG"
echo "📝 Log: $LOG_FILE"

# Verifica se já está rodando
if [ -f "$PID_FILE" ]; then
    OLD_PID=$(cat "$PID_FILE")
    if kill -0 "$OLD_PID" 2>/dev/null; then
        echo "❌ Pipeline já está rodando (PID: $OLD_PID)"
        echo "   Para parar: kill $OLD_PID"
        echo "   Para monitorar: tail -f $LOG_FILE"
        exit 1
    else
        echo "🧹 Removendo PID file órfão"
        rm -f "$PID_FILE"
    fi
fi

# Verifica se config existe
if [ ! -f "$CONFIG" ]; then
    echo "❌ Arquivo de configuração não encontrado: $CONFIG"
    exit 1
fi

# Executa em background com nohup
echo "🚀 Iniciando pipeline em background..."
echo "📝 Logs serão salvos em: $LOG_FILE"
echo "🔍 Para monitorar: tail -f $LOG_FILE"
echo "⏹️  Para parar: kill \$(cat $PID_FILE)"

# Executa o pipeline
nohup ./genomes_analyzer.py --config "$CONFIG" > "$LOG_FILE" 2>&1 &
PIPELINE_PID=$!

# Salva PID
echo "$PIPELINE_PID" > "$PID_FILE"

echo "✅ Pipeline iniciado!"
echo "📊 PID: $PIPELINE_PID"
echo "📝 Log: $LOG_FILE"
echo ""
echo "📋 Comandos úteis:"
echo "   Monitorar:     tail -f $LOG_FILE"
echo "   Monitorar live: watch 'tail -20 $LOG_FILE'"
echo "   Status:        ps -p $PIPELINE_PID"
echo "   Parar:         kill $PIPELINE_PID"
echo "   Recursos:      htop -p $PIPELINE_PID"
echo ""
echo "🎯 Configuração ativada!"
echo ""

# Mostra início do log
echo "📖 Primeiras linhas do log:"
sleep 2
head -20 "$LOG_FILE" || echo "Log ainda não iniciado..."
