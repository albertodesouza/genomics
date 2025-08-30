#!/bin/bash
# setup_monster_256gb.sh - Configurações do sistema para máquina monster 256GB

echo "🚀 Configurando sistema para máquina monster (128 cores + 256GB RAM)"

# Verifica RAM disponível
TOTAL_RAM=$(free -g | awk 'NR==2{print $2}')
echo "💾 RAM total detectada: ${TOTAL_RAM}GB"

if [ "$TOTAL_RAM" -lt 200 ]; then
    echo "⚠️  Aviso: RAM detectada ($TOTAL_RAM GB) é menor que 256GB"
    echo "   Configuração pode precisar de ajustes"
fi

# Configura /dev/shm para usar mais RAM (padrão é 50%)
echo "🔧 Configurando /dev/shm para usar até 128GB..."
sudo mount -o remount,size=128G /dev/shm
echo "✅ /dev/shm configurado para 128GB"

# Verifica espaço em /dev/shm
SHM_SIZE=$(df -h /dev/shm | awk 'NR==2{print $2}')
echo "💿 /dev/shm disponível: $SHM_SIZE"

# Otimizações do kernel para I/O intensivo
echo "⚙️  Aplicando otimizações do kernel..."

# Aumenta limites de file descriptors
echo "fs.file-max = 2000000" | sudo tee -a /etc/sysctl.conf
echo "* soft nofile 1000000" | sudo tee -a /etc/security/limits.conf
echo "* hard nofile 1000000" | sudo tee -a /etc/security/limits.conf

# Otimizações para I/O
echo "vm.dirty_ratio = 15" | sudo tee -a /etc/sysctl.conf
echo "vm.dirty_background_ratio = 5" | sudo tee -a /etc/sysctl.conf
echo "vm.swappiness = 1" | sudo tee -a /etc/sysctl.conf

# Aplica configurações
sudo sysctl -p

echo ""
echo "🎯 Configurações aplicadas para máquina monster:"
echo "   • RAM disk /dev/shm: 128GB"
echo "   • File descriptors: 1M"
echo "   • I/O otimizado para workloads intensivos"
echo "   • Swappiness mínimo (RAM abundante)"
echo ""
echo "💡 Para reverter configurações:"
echo "   sudo mount -o remount,size=50% /dev/shm"
echo ""
echo "✅ Sistema otimizado para pipeline genomics!"
echo "🚀 Agora execute: ./run_monster_background.sh"
