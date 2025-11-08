# 🚀 Guia de Instalação - Máquina Monster (128 cores + 256GB RAM)

## 📋 Pré-requisitos

- **SO**: Linux (Ubuntu/CentOS/RHEL)
- **Hardware**: 128+ cores, 256GB+ RAM
- **Espaço**: 5TB+ livre
- **Rede**: Conexão estável para downloads

## 🔧 Instalação Passo a Passo

### **1. Transferir Arquivos para Monster**

```bash
# Na sua máquina local, compacte os arquivos:
tar -czf genomics_pipeline.tar.gz *.py *.sh *.yaml *.md

# Transfira para a monster:
scp genomics_pipeline.tar.gz usuario@monster:/home/usuario/

# Na monster, extraia:
ssh usuario@monster
cd /home/usuario
tar -xzf genomics_pipeline.tar.gz
cd genomics
```

### **2. Instalar Conda/Mamba**

```bash
# Torna script executável
chmod +x scripts/install_conda_universal.sh

# Instala conda/mamba automaticamente
./scripts/install_conda_universal.sh

# Reinicia terminal ou carrega configurações
source ~/.bashrc
```

### **3. Instalar Ambiente Genomics**

```bash
# Torna scripts executáveis
chmod +x scripts/*.sh

# Instala ambiente genomics
./scripts/install_genomics_env.sh

# Instala VEP
source scripts/vep_install.sh
```

### **4. Configurar Sistema Monster**

```bash
# Otimiza sistema para 256GB RAM
./scripts/setup_monster_256gb.sh
```

### **5. Testar Instalação**

```bash
# Testa ambiente
source scripts/start_genomics_universal.sh

# Verifica ferramentas
bwa-mem2 version
samtools --version
bcftools --version
vep --help | head -10
```

## 🚀 Execução

### **Execução em Background:**
```bash
./scripts/run_monster_background.sh
```

### **Monitoramento:**
```bash
./scripts/monitor_monster.sh
```

## 🔧 Troubleshooting

### **Se conda não for encontrado:**
```bash
# Verifica instalação
ls -la ~/miniforge3/
which conda

# Recarrega configurações
source ~/.bashrc
```

### **Se ambiente genomics não existir:**
```bash
# Lista ambientes
conda env list

# Recria ambiente
./scripts/install_genomics_env.sh
```

### **Se VEP não funcionar:**
```bash
# Reinstala VEP
source scripts/vep_install.sh

# Verifica cache
ls -la /dados/vep_cache/
```

## 📊 Performance Esperada

**Máquina Monster vs Normal:**

| Etapa | Normal (16 cores) | Monster (128 cores) | Speedup |
|-------|-------------------|---------------------|---------|
| **Download** | 2-4h | 30-60min | **4-6x** |
| **Alinhamento** | 6-8h | 45-90min | **6-8x** |
| **Variant Calling** | 4-6h | 30-45min | **8-10x** |
| **VEP** | 8-12h | 60-90min | **6-8x** |
| **Total** | **20-30h** | **3-4h** | **6-8x** |

## 💾 Uso de Recursos

**RAM (256GB):**
- Alinhamento: ~150GB
- GATK: ~128GB heap
- VEP: ~50GB cache
- RAM disk: ~128GB

**CPU (128 cores):**
- 96 cores alinhamento
- 32 cores variant calling
- 96 cores VEP

**Disco:**
- Base: ~3TB (dados 30× completos)
- Temp: ~1TB (RAM disk)
- Cache: ~30GB (VEP)

## 🎯 Configurações Específicas

**Para máquinas com diferentes specs:**

### **64 cores + 128GB RAM:**
```yaml
aln_threads: 48
bcf_max_parallel: 16
vep_fork: 48
hc_java_mem_gb: 64
```

### **32 cores + 64GB RAM:**
```yaml
aln_threads: 24
bcf_max_parallel: 8
vep_fork: 24
hc_java_mem_gb: 32
```

## 🆘 Suporte

**Se encontrar problemas:**

1. **Verifique logs**: `tail -f pipeline_monster_*.log`
2. **Monitore recursos**: `htop`
3. **Verifique espaço**: `df -h`
4. **Teste ferramentas**: `source start_genomics_universal.sh`
