# 🚀 Quick Start - Projeto de Longevidade

Guia rápido para começar com o projeto de descoberta de marcadores de longevidade.

---

## ⚡ Setup Rápido (5 minutos)

```bash
# 1. Ativar ambiente
conda activate genomics

# 2. Instalar dependências adicionais
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install scikit-learn

# 3. Configurar API key do AlphaGenome
nano configs/default.yaml
# Editar linha: api_key: "YOUR_ALPHAGENOME_API_KEY"
```

---

## 🔄 Pipeline Completo em 3 Comandos

### 1. Construir Dataset

```bash
python neural_longevity_dataset.py --config configs/default.yaml
```

**O que acontece**:
- Download de lista de amostras do 1000 Genomes
- Seleção de 30 amostras "longevas" (posições 0-30)
- Seleção de 30 amostras "não-longevas" (posições 500-530)
- Identificação de 10 variantes da primeira amostra longeva
- Extração de sequências FASTA (2048 bp) centradas nas variantes
- Processamento com AlphaGenome
- Construção de dataset PyTorch (train/val/test splits)

**Tempo estimado**: 2-4 horas (depende da API AlphaGenome)

**Saídas**:
```
longevity_dataset/
├── sequences/              # Arquivos FASTA
├── alphagenome_cache/      # Cache de predições
├── central_points.json     # Pontos centrais
├── train.pkl              # Dataset treino
├── val.pkl                # Dataset validação
└── test.pkl               # Dataset teste
```

### 2. Treinar Modelo

```bash
python longevity_train.py --config longevity_train_config.yaml
```

**O que acontece**:
- Carrega datasets PyTorch
- Cria modelo CNN
- Treina por 100 épocas (com early stopping)
- Valida a cada época
- Salva melhor modelo

**Tempo estimado**: 30-60 minutos (CPU), 5-10 minutos (GPU)

**Saídas**:
```
longevity_models/
├── best_model.pt          # Melhor modelo
├── history.json           # Histórico de treino
└── test_results.json      # Resultados no teste
```

### 3. Avaliar Resultados

```bash
python longevity_train.py \
  --config longevity_train_config.yaml \
  --test-only
```

**O que acontece**:
- Carrega melhor modelo
- Avalia no conjunto de teste
- Exibe métricas: Accuracy, Precision, Recall, F1, AUC-ROC
- Matriz de confusão

---

## 📊 Exemplo de Uso do Dataset

```python
import torch
from torch.utils.data import DataLoader
from neural_longevity_dataset import LongevityDataset

# Carregar dataset
train_dataset = LongevityDataset('longevity_dataset/train.pkl')

# Informações
print(f"Total de amostras: {len(train_dataset)}")
print(f"Distribuição de classes: {train_dataset.get_class_distribution()}")

# Inspecionar uma amostra
features, label = train_dataset[0]
print(f"\nSequence shape: {features['sequence'].shape}")      # (4, 2048)
print(f"Position: {features['position']}")                    # scalar
print(f"AlphaGenome features: {features['alphagenome'].shape}") # (35,)
print(f"Label: {label}")  # 0 ou 1
print(f"Metadata: {features['metadata']}")

# Criar DataLoader
loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

for batch_features, batch_labels in loader:
    print(f"Batch de {len(batch_labels)} amostras")
    break
```

---

## 🧪 Modo Desenvolvimento/Teste

Para testar rapidamente sem processar tudo:

### 1. Dry-run (Simulação)

```bash
python neural_longevity_dataset.py \
  --config configs/default.yaml \
  --dry-run
```

Mostra o que seria feito sem executar.

### 2. Limitar Amostras

Edite `configs/default.yaml`:

```yaml
debug:
  limit_samples: 5  # Apenas 5 amostras de cada grupo
```

### 3. Usar Dados Simulados

Se não tiver VCFs reais, o pipeline cria pontos centrais simulados automaticamente.

---

## 🔧 Configurações Importantes

### configs/default.yaml

```yaml
# Tamanho da janela de DNA (deve ser suportado pelo AlphaGenome)
sequence_extraction:
  window_size: 2048  # Opções: 2048, 16384, 131072, 524288, 1048576

# Número de pontos centrais (variantes)
variant_selection:
  n_central_points: 10  # Comece com 10, aumente depois

# Outputs do AlphaGenome
alphagenome:
  outputs:
    - "RNA_SEQ"
    - "CAGE"
    - "ATAC"
    - "CHIP_HISTONE"
    - "CHIP_TF"
    - "DNASE"
    - "PROCAP"
```

### longevity_train_config.yaml

```yaml
training:
  epochs: 100
  batch_size: 32
  learning_rate: 0.001
```

---

## 📈 Monitorar Treinamento

Durante o treino, você verá:

```
Época 1/100
Treinando: 100%|███████████| 30/30 [00:15<00:00]
Validando: 100%|███████████| 10/10 [00:03<00:00]

┏━━━━━━━━━━━━━┳━━━━━━━━┓
┃ Métrica     ┃ Valor  ┃
┡━━━━━━━━━━━━━╇━━━━━━━━┩
│ Train Loss  │ 0.6845 │
│ Val Loss    │ 0.6521 │
│ Accuracy    │ 0.6333 │
│ Precision   │ 0.6500 │
│ Recall      │ 0.6000 │
│ F1-Score    │ 0.6240 │
│ AUC-ROC     │ 0.6880 │
└─────────────┴────────┘
✓ Melhor modelo salvo!
```

---

## 🎯 Resultados Esperados

### Baseline (Dados Simulados)

Com dados simulados (sem correlação real com longevidade):
- **Accuracy**: ~50-60% (aleatório)
- **AUC-ROC**: ~0.5-0.6

### Com Dados Reais (Objetivo)

Com dados reais de pessoas longevas:
- **Accuracy**: >70% (objetivo)
- **AUC-ROC**: >0.75 (objetivo)
- **F1-Score**: >0.70

---

## 🔍 Troubleshooting

### Erro: API key não configurada

```
✗ API key do AlphaGenome não configurada!
```

**Solução**: Edite `configs/default.yaml` e adicione sua API key.

### Erro: VCF não encontrado

```
⚠ VCF não encontrado: vcfs_longevous/SAMPLE_0000.vcf.gz
```

**Solução**: O pipeline criará pontos simulados automaticamente. Para usar dados reais:
1. Processe genomas com `genomes_analyzer.py`
2. Coloque VCFs em `longevity_dataset/vcfs_longevous/` e `vcfs_non_longevous/`

### Erro: Memória insuficiente

```
CUDA out of memory
```

**Solução**: Reduza batch_size em `longevity_train_config.yaml`:
```yaml
training:
  batch_size: 16  # Era 32
```

### Processamento muito lento

**Solução**: Use cache do AlphaGenome (já ativado por padrão):
```yaml
alphagenome:
  cache_results: true
```

---

## 📚 Próximos Passos

Após validar o pipeline com dados simulados:

1. **Obter Dados Reais de Longevidade**
   - Bases de dados públicas
   - Colaborações com instituições
   
2. **Refinar Pontos Centrais**
   - Análise de importância de features
   - Seleção dos top-K pontos mais discriminativos
   
3. **Experimentos com Arquitetura**
   - Testar window_size diferentes
   - Adicionar mais outputs do AlphaGenome
   - Experimentar com Transformers
   
4. **Validação Biológica**
   - Verificar se marcadores têm sentido biológico
   - Consultar literatura sobre longevidade
   - Colaborar com biólogos/médicos

---

## 🆘 Ajuda

- **Documentação completa**: [docs/PROJECT.md](docs/PROJECT.md)
- **Neural Module**: [../neural_module/README.md](../neural_module/README.md)
- **Pipeline genômico**: [../README.md](../README.md)

---

**Bom trabalho! 🚀🧬**

