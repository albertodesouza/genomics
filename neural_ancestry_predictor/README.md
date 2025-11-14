# Neural Ancestry Predictor

> **🧬 Predição de Ancestralidade usando Redes Neurais e Dados AlphaGenome**

Este módulo implementa uma rede neural configurável via YAML que prediz ancestralidade (superpopulation, population ou FROG likelihood) a partir de predições AlphaGenome armazenadas em um dataset PyTorch.

## 📑 Índice

- [Visão Geral](#visão-geral)
- [Instalação](#instalação)
- [Uso Rápido](#uso-rápido)
- [Configuração](#configuração)
- [Arquitetura](#arquitetura)
- [Processamento de Dados](#processamento-de-dados)
- [Treinamento](#treinamento)
- [Teste e Avaliação](#teste-e-avaliação)
- [Weights & Biases](#weights--biases)
- [Ajuste de Hiperparâmetros](#ajuste-de-hiperparâmetros)
- [FAQ](#faq)

---

## Visão Geral

### O que este módulo faz?

O **Neural Ancestry Predictor** treina uma rede neural para prever a ancestralidade genética de indivíduos usando:

- **Entrada**: Predições AlphaGenome (ex: ATAC-seq, RNA-seq) de janelas genômicas
- **Saída**: Superpopulação (AFR, AMR, EAS, EUR, SAS), População (26 classes) ou FROG likelihood (150 valores)

### Características

- ✅ **Totalmente configurável via YAML**
- ✅ **Suporta múltiplos targets** (superpopulation, population, FROG likelihood)
- ✅ **Processamento flexível** de janelas, haplótipos e outputs AlphaGenome
- ✅ **Integração com Weights & Biases** para tracking e visualização
- ✅ **Checkpointing automático** para salvar modelos treinados
- ✅ **Métricas detalhadas** (accuracy, precision, recall, F1, confusion matrix)
- ✅ **Normalização automática** com cache de parâmetros

---

## Instalação

### 1. Dependências

```bash
# Navegar para o diretório
cd genomics/neural_ancestry_predictor

# Instalar dependências Python
pip install torch numpy pandas pyyaml scikit-learn rich

# Opcional: Weights & Biases (para tracking)
pip install wandb
wandb login  # Autenticar com sua conta W&B
```

### 2. Dataset

Este módulo requer um dataset PyTorch criado pelo `build_non_longevous_dataset`:

```bash
# Exemplo: dataset em /dados/GENOMICS_DATA/top3/non_longevous_results
# Deve conter:
#   - individuals/
#   - dataset_metadata.json
```

Consulte `build_non_longevous_dataset/docs/PYTORCH_DATASET.md` para mais informações.

---

## Uso Rápido

### Treinar Modelo

```bash
cd neural_ancestry_predictor
python3 neural_ancestry_predictor.py --config configs/default.yaml
```

### Testar Modelo

```bash
python3 neural_ancestry_predictor.py --config configs/default.yaml --mode test
```

### Exemplo de Saída

```
╭─────────────────────────────────────────╮
│ 🧬 Genomics                             │
│                                         │
│ Neural Ancestry Predictor               │
│ Modo: train                             │
│ Target: superpopulation                 │
│ Config: configs/default.yaml            │
╰─────────────────────────────────────────╯

Device: cuda
[INFO] GenomicLongevityDataset inicializado:
  • Dataset: non_longevous_1000g
  • Indivíduos: 78
  • Load predictions: True
  • Load sequences: False

Computando parâmetros de normalização...
Normalização: mean=0.123456, std=0.654321

Dataset split:
  • Treino: 54 amostras
  • Validação: 12 amostras
  • Teste: 12 amostras

Modelo criado:
  • Input size: 11000
  • Hidden layers: [128, 64]
  • Output size: 5
  • Activation: relu
  • Dropout: 0.2
  • Total parameters: 1,415,237

╭─────────────────────────────────────────╮
│ Iniciando Treinamento                   │
│                                         │
│ Épocas: 100                             │
│ Batch size: 16                          │
│ Learning rate: 0.001                    │
╰─────────────────────────────────────────╯

Validação - Época 5: Loss=0.8234, Accuracy=0.7500
✓ Checkpoint salvo: models/best_accuracy.pt

...

✓ Treinamento concluído!
```

---

## Configuração

### Estrutura do Arquivo YAML

O arquivo `configs/default.yaml` contém **8 seções principais**:

#### A) Dataset Input Parameters

Controla **o que** carregar do dataset e **como** processar:

```yaml
dataset_input:
  dataset_dir: "/path/to/dataset"           # Caminho do dataset PyTorch
  alphagenome_outputs: ["ATAC"]             # Quais outputs usar (RNA_SEQ, ATAC, CAGE, etc.)
  haplotype_mode: "H1+H2"                   # "H1", "H2" ou "H1+H2"
  window_center_size: 100                   # Tamanho do trecho central (bases)
  downsample_factor: 1                      # Fator de downsampling (1 = sem)
```

**Impacto na Dimensionalidade:**

- Cada janela tem ~1M bases por output
- `window_center_size=100` → extrai 100 bases do centro
- `downsample_factor=2` → usa 1 a cada 2 bases (reduz para 50)
- `haplotype_mode="H1+H2"` → dobra o tamanho (2 haplótipos)
- Dimensão final = `n_windows × n_outputs × n_haplotypes × (window_center_size / downsample_factor)`

**Exemplo:**
- 55 janelas (SNPs)
- 1 output (ATAC)
- 2 haplótipos (H1+H2)
- window_center_size=100, downsample_factor=1
- **Dimensão = 55 × 1 × 2 × 100 = 11,000 features**

#### B) Output Parameters

Define **o que** a rede deve prever:

```yaml
output:
  prediction_target: "superpopulation"  # "superpopulation", "population" ou "frog_likelihood"
```

| Target | Tipo | Classes | Dificuldade |
|--------|------|---------|-------------|
| `superpopulation` | Classificação | 5 (AFR, AMR, EAS, EUR, SAS) | Fácil ⭐ |
| `population` | Classificação | 26 | Média ⭐⭐ |
| `frog_likelihood` | Regressão | 150 valores | Difícil ⭐⭐⭐ |

#### C) Model Architecture

Define a **arquitetura** da rede:

```yaml
model:
  hidden_layers: [128, 64]    # Lista de neurônios por camada oculta
  activation: "relu"          # "relu", "tanh" ou "sigmoid"
  dropout_rate: 0.2           # Taxa de dropout (0.0 a 1.0)
```

**Arquitetura resultante:**

```
Input (11000) → Dense(128) → ReLU → Dropout(0.2) →
Dense(64) → ReLU → Dropout(0.2) →
Dense(5) → Softmax
```

#### D) Training Parameters

Controla o **processo de treinamento**:

```yaml
training:
  optimizer: "adam"              # "adam", "adamw" ou "sgd"
  learning_rate: 0.001           # Taxa de aprendizado
  loss_function: "cross_entropy" # "cross_entropy" ou "mse"
  batch_size: 16                 # Número de amostras por batch
  num_epochs: 100                # Número de épocas
  validation_frequency: 5        # Validar a cada N épocas
```

#### E) Data Split

Define **divisão do dataset**:

```yaml
data_split:
  train_split: 0.7      # 70% para treino
  val_split: 0.15       # 15% para validação
  test_split: 0.15      # 15% para teste
  random_seed: 42       # Seed para reprodutibilidade
```

#### F) Weights & Biases

Configuração de **tracking e visualização**:

```yaml
wandb:
  use_wandb: false                          # Habilitar W&B
  project_name: "neural-ancestry-predictor" # Nome do projeto
  run_name: null                            # Nome do run (auto se null)
  log_frequency: 10                         # Log a cada N batches
```

#### G) Checkpointing

Controla **salvamento de modelos**:

```yaml
checkpointing:
  checkpoint_dir: "models"       # Diretório para checkpoints
  save_frequency: 10             # Salvar a cada N épocas
  load_checkpoint: null          # Caminho para checkpoint existente
```

#### H) Mode

Define **modo de operação**:

```yaml
mode: "train"  # "train" ou "test"
```

---

## Arquitetura

### Visão Geral

```
┌─────────────────────────────────────────────────────────────┐
│                    NEURAL ANCESTRY PREDICTOR                 │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Input: AlphaGenome Predictions                             │
│  ┌──────────────────────────────────────────┐               │
│  │ Window 1 (SNP/Gene)                      │               │
│  │  ├─ H1: [ATAC: 100 bases]                │               │
│  │  └─ H2: [ATAC: 100 bases]                │               │
│  │ Window 2                                  │               │
│  │  ├─ H1: [ATAC: 100 bases]                │               │
│  │  └─ H2: [ATAC: 100 bases]                │               │
│  │ ...                                       │               │
│  │ Window 55                                 │               │
│  │  ├─ H1: [ATAC: 100 bases]                │               │
│  │  └─ H2: [ATAC: 100 bases]                │               │
│  └──────────────────────────────────────────┘               │
│           ↓ Concatenação + Normalização                     │
│  ┌──────────────────────────────────────────┐               │
│  │ Feature Vector [11000 elementos]         │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  ┌──────────────────────────────────────────┐               │
│  │ Dense Layer (128 neurons)                │               │
│  │ ReLU Activation                          │               │
│  │ Dropout (0.2)                            │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  ┌──────────────────────────────────────────┐               │
│  │ Dense Layer (64 neurons)                 │               │
│  │ ReLU Activation                          │               │
│  │ Dropout (0.2)                            │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  ┌──────────────────────────────────────────┐               │
│  │ Output Layer (5 neurons)                 │               │
│  │ Softmax                                  │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  Output: [AFR, AMR, EAS, EUR, SAS] probabilities           │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Componentes

1. **ProcessedGenomicDataset**: Dataset wrapper que:
   - Carrega dados do `GenomicLongevityDataset`
   - Extrai trecho central das janelas
   - Aplica downsampling
   - Combina haplótipos
   - Normaliza (z-score)
   - Cache de parâmetros de normalização

2. **AncestryPredictor**: Modelo PyTorch que:
   - Construção dinâmica baseada em config
   - Camadas totalmente conectadas (Dense)
   - Dropout para regularização
   - Softmax para classificação ou linear para regressão

3. **Trainer**: Gerencia treinamento:
   - Loop de treino com progress bars
   - Validação periódica
   - Checkpointing automático
   - Logging no W&B

4. **Tester**: Gerencia teste:
   - Inferência no conjunto de teste
   - Métricas detalhadas
   - Confusion matrix
   - Classification report

---

## Processamento de Dados

### Pipeline de Processamento

```
Dataset Original → Extração → Downsampling → Combinação → Normalização → Tensor
                   de Centro                   Haplótipos
```

#### 1. Extração de Trecho Central

Cada janela tem ~1M bases. Extraímos o trecho central:

```python
# window_center_size = 100
# Array original: [1000000 elementos]
center_idx = 500000
start = center_idx - 50  # 499950
end = center_idx + 50    # 500050
extracted = array[start:end]  # [100 elementos]
```

**Por que centro?** Assume que região central é mais relevante (próxima ao gene/SNP).

#### 2. Downsampling

Reduz ainda mais a dimensionalidade:

```python
# downsample_factor = 2
downsampled = extracted[::2]  # [50 elementos]
```

#### 3. Combinação de Haplótipos

```python
# haplotype_mode = "H1+H2"
features = concatenate([H1_features, H2_features])

# haplotype_mode = "H1"
features = H1_features
```

#### 4. Normalização

Normalização z-score usando toda a base:

```python
# Computado uma vez no início
mean = mean(all_training_data)
std = std(all_training_data)

# Aplicado a cada amostra
normalized = (features - mean) / std
```

Parâmetros salvos em `models/normalization_params.json` para reutilização.

---

## Treinamento

### Executar Treinamento

```bash
python3 neural_ancestry_predictor.py --config configs/default.yaml
```

### Durante o Treinamento

O programa irá:

1. **Carregar dataset** e imprimir sumário
2. **Computar normalização** (pode demorar alguns minutos)
3. **Dividir dados** em treino/validação/teste
4. **Criar modelo** e imprimir arquitetura
5. **Treinar** com progress bars por época
6. **Validar** a cada N épocas
7. **Salvar checkpoints**:
   - `best_loss.pt`: Melhor loss de validação
   - `best_accuracy.pt`: Melhor accuracy de validação
   - `epoch_N.pt`: Checkpoints periódicos

### Monitoramento

**Terminal:**
- Progress bars por época
- Loss e accuracy de validação
- Avisos e erros

**Arquivos:**
- `models/training_history.json`: Histórico completo
- `models/normalization_params.json`: Parâmetros de normalização

**Weights & Biases** (se habilitado):
- Gráficos de loss em tempo real
- Métricas de validação
- Histogramas de gradientes
- Comparação entre runs

### Continuar Treinamento

Para continuar de um checkpoint:

```yaml
# configs/default.yaml
checkpointing:
  load_checkpoint: "models/epoch_50.pt"
```

---

## Teste e Avaliação

### Executar Teste

```bash
python3 neural_ancestry_predictor.py --config configs/default.yaml --mode test
```

Ou configure no YAML:

```yaml
mode: "test"
checkpointing:
  load_checkpoint: "models/best_accuracy.pt"
```

### Resultados

O teste gera:

**1. Métricas Gerais:**

```
╔═══════════════════════════════════════╗
║        Métricas de Performance        ║
╠═══════════════════════════════════════╣
║ Accuracy           │ 0.9167           ║
║ Precision (weighted)│ 0.9250          ║
║ Recall (weighted)  │ 0.9167           ║
║ F1-Score (weighted)│ 0.9183           ║
╚═══════════════════════════════════════╝
```

**2. Classification Report:**

```
              precision    recall  f1-score   support

         AFR       0.92      0.95      0.93        20
         AMR       0.88      0.85      0.86        13
         EAS       0.95      0.93      0.94        15
         EUR       0.90      0.93      0.91        15
         SAS       0.93      0.90      0.91        15

    accuracy                           0.92        78
   macro avg       0.92      0.91      0.91        78
weighted avg       0.92      0.92      0.92        78
```

**3. Confusion Matrix:**

```
╔════════════════════════════════════════════════╗
║              Confusion Matrix                  ║
╠════════════════════════════════════════════════╣
║ True \ Pred │  AFR  │  AMR  │  EAS  │  EUR  │  SAS  ║
║ AFR         │   19  │    1  │    0  │    0  │    0  ║
║ AMR         │    1  │   11  │    0  │    1  │    0  ║
║ EAS         │    0  │    0  │   14  │    1  │    0  ║
║ EUR         │    0  │    1  │    0  │   14  │    0  ║
║ SAS         │    0  │    0  │    1  │    0  │   14  ║
╚════════════════════════════════════════════════╝
```

### Interpretação

- **Accuracy**: % de predições corretas
- **Precision**: % de predições positivas corretas
- **Recall**: % de casos positivos identificados
- **F1-Score**: Média harmônica de precision e recall
- **Confusion Matrix**: Onde o modelo erra

**Exemplo de análise:**
- AFR: Alta recall (0.95) → identifica bem africanos
- AMR: Menor precision (0.88) → às vezes confunde com outras
- Diagonal forte → modelo bem calibrado

---

## Weights & Biases

### Configurar W&B

```bash
# Instalar
pip install wandb

# Autenticar
wandb login

# Habilitar no config
```

```yaml
wandb:
  use_wandb: true
  project_name: "neural-ancestry-predictor"
  run_name: "experiment-atac-h1h2-100bases"  # Opcional
```

### Visualizações Disponíveis

1. **Loss Curves**: Train vs Validation loss
2. **Accuracy**: Evolução da accuracy
3. **Confusion Matrix**: Matriz interativa
4. **Gradients**: Histogramas de gradientes
5. **Parameters**: Distribuição de pesos
6. **System Metrics**: GPU, CPU, RAM

### Comparar Experimentos

No dashboard do W&B, você pode:
- Sobrepor gráficos de múltiplos runs
- Filtrar por hiperparâmetros
- Gerar tabelas de comparação
- Exportar gráficos para papers (PNG, SVG, PDF)

---

## Ajuste de Hiperparâmetros

### Dimensionalidade muito alta?

**Problema**: Treino muito lento, memória insuficiente

**Soluções**:
```yaml
dataset_input:
  window_center_size: 50        # Reduzir de 100 para 50
  downsample_factor: 2          # Usar 1 a cada 2 bases
  haplotype_mode: "H1"          # Usar apenas um haplótipo
  alphagenome_outputs: ["ATAC"] # Usar apenas 1 output
```

### Underfitting (loss alto)?

**Problema**: Modelo não consegue aprender padrões

**Soluções**:
```yaml
model:
  hidden_layers: [256, 128, 64]  # Mais camadas/neurônios
  dropout_rate: 0.0               # Remover regularização

training:
  num_epochs: 200                 # Treinar por mais tempo
```

### Overfitting (val_loss > train_loss)?

**Problema**: Modelo memoriza treino mas não generaliza

**Soluções**:
```yaml
model:
  dropout_rate: 0.5               # Aumentar dropout
  hidden_layers: [64, 32]         # Reduzir capacidade

training:
  learning_rate: 0.0001           # Learning rate menor
```

### Treino instável (loss oscila)?

**Problema**: Gradientes instáveis

**Soluções**:
```yaml
training:
  learning_rate: 0.0001           # LR menor
  optimizer: "adamw"              # Tentar outro otimizador
  batch_size: 32                  # Batch maior
```

### Convergência lenta?

**Problema**: Treino demora muito

**Soluções**:
```yaml
training:
  learning_rate: 0.01             # LR maior (cuidado!)
  batch_size: 64                  # Batch maior
  optimizer: "adam"               # Adam geralmente mais rápido que SGD
```

---

## FAQ

### Q: Quanto tempo leva o treinamento?

**A**: Depende de:
- **Dataset size**: 78 amostras → minutos; 1000 amostras → horas
- **Dimensionalidade**: 100 bases → rápido; 10000 bases → lento
- **Hardware**: GPU → 10-100x mais rápido que CPU
- **Épocas**: 100 épocas → proporcional

**Estimativas** (78 amostras, 100 bases, GPU):
- Normalização: ~30 segundos
- Época: ~5 segundos
- 100 épocas: ~8 minutos

### Q: Qual target devo usar?

**A**: Recomendações:
- **Iniciante**: `superpopulation` (5 classes, mais fácil)
- **Intermediário**: `population` (26 classes)
- **Avançado**: `frog_likelihood` (regressão, 150 valores)

### Q: Preciso de GPU?

**A**: Não é obrigatório, mas **altamente recomendado**:
- CPU: Funciona, mas ~50-100x mais lento
- GPU: Nvidia com CUDA (RTX 3060 ou superior ideal)

Para instalar PyTorch com GPU:
```bash
# CUDA 11.8
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

### Q: Como interpretar a accuracy?

**A**: Depende do baseline:
- **Random guessing** (5 classes): 20%
- **Bom modelo**: 70-80%
- **Excelente modelo**: 85-95%
- **Perfeito**: 100% (cuidado com overfitting!)

Compare sempre com validação E teste.

### Q: Posso usar múltiplos outputs AlphaGenome?

**A**: Sim! Aumenta dimensionalidade mas pode melhorar desempenho:

```yaml
dataset_input:
  alphagenome_outputs:
    - "ATAC"
    - "RNA_SEQ"
    - "CAGE"
```

Dimensão cresce linearmente: 1 output → 11k features; 3 outputs → 33k features.

### Q: Como adicionar mais genes/SNPs ao dataset?

**A**: Recrie o dataset com `build_non_longevous_dataset`:

```yaml
# build_non_longevous_dataset/configs/custom.yaml
build_window_params:
  mode: "snp"
  snp:
    snp_list_file: "path/to/your_snps.txt"  # Adicionar mais SNPs
```

Mais janelas → maior dimensionalidade → pode melhorar ou piorar (curse of dimensionality).

### Q: Erro "CUDA out of memory"?

**A**: Reduza uso de memória:

```yaml
training:
  batch_size: 4          # Reduzir batch
  
dataset_input:
  window_center_size: 50  # Reduzir dimensão
  downsample_factor: 2    # Aumentar downsampling
```

Ou use CPU:
```bash
# Forçar CPU
export CUDA_VISIBLE_DEVICES=""
python3 neural_ancestry_predictor.py --config configs/default.yaml
```

### Q: Como exportar gráficos para paper?

**A**: 

**Opção 1: Weights & Biases**
- No dashboard, clicar em gráfico → "Export" → PNG/SVG
- Alta qualidade, ideal para publicação

**Opção 2: Programaticamente**
```python
import matplotlib.pyplot as plt
import json

# Carregar histórico
with open('models/training_history.json') as f:
    history = json.load(f)

# Plotar
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(history['epoch'], history['train_loss'], label='Train')
plt.plot(history['epoch'], history['val_loss'], label='Validation')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.savefig('loss_curve.png', dpi=300, bbox_inches='tight')
```

---

## Estrutura de Arquivos

```
neural_ancestry_predictor/
├── neural_ancestry_predictor.py    # Programa principal
├── configs/
│   └── default.yaml                 # Configuração padrão
├── models/                          # Checkpoints (criado automaticamente)
│   ├── best_loss.pt
│   ├── best_accuracy.pt
│   ├── epoch_10.pt
│   ├── epoch_20.pt
│   ├── normalization_params.json
│   └── training_history.json
└── README.md                        # Esta documentação
```

---

## Referências

- **PyTorch**: https://pytorch.org/
- **Weights & Biases**: https://wandb.ai/
- **AlphaGenome**: https://alphagenome.ai/
- **1000 Genomes**: http://www.internationalgenome.org/
- **FROGAncestryCalc**: Ancestry inference via AISNPs

---

## Suporte

Para problemas ou questões:

1. Verifique este README
2. Consulte `build_non_longevous_dataset/docs/PYTORCH_DATASET.md`
3. Execute com modo debug: adicione prints no código
4. Verifique logs do W&B (se habilitado)

**Autor**: Alberto F. De Souza (via ChatGPT)  
**Data**: 2025-11-14  
**Versão**: 1.0

---

## Changelog

### v1.0 (2025-11-14)
- ✨ Implementação inicial
- ✨ Suporte para superpopulation, population e FROG likelihood
- ✨ Integração com Weights & Biases
- ✨ Processamento configurável de janelas e haplótipos
- ✨ Checkpointing e normalização automática
- ✨ Métricas detalhadas e visualizações

