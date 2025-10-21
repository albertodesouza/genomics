# 🧬 IA Neuro-Simbólica para Descoberta de Marcadores Genéticos e Epigenéticos de Longevidade

## 📖 Visão Geral do Projeto

Este projeto implementa um pipeline completo para descoberta de marcadores de longevidade usando:
- **Análise Genômica Tradicional** (`genomes_analyzer.py`)
- **IA Generativa** (AlphaGenome da Google DeepMind)
- **Aprendizado Profundo** (PyTorch para classificação)

### 🎯 Objetivo

Treinar uma rede neural profunda capaz de distinguir, a partir do DNA, se uma pessoa será ou foi longeva, usando como features as predições funcionais do AlphaGenome sobre regiões genômicas contendo variantes.

---

## 🔄 Pipeline Completo

```
┌─────────────────────────────────────────────────────────────┐
│ 1. COLETA DE DADOS                                          │
├─────────────────────────────────────────────────────────────┤
│ • Genomas de pessoas longevas                               │
│ • Genomas de pessoas não-longevas (controle)               │
│ • Fonte inicial: 1000 Genomes Project                       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. IDENTIFICAÇÃO DE VARIANTES                               │
├─────────────────────────────────────────────────────────────┤
│ • Comparação com genoma de referência (GRCh38)              │
│ • Extração de SNVs, INDELs                                  │
│ • Filtros de qualidade                                      │
│ • Seleção de pontos centrais                                │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. EXTRAÇÃO DE SEQUÊNCIAS                                   │
├─────────────────────────────────────────────────────────────┤
│ • Janelas de DNA centradas nas variantes                    │
│ • Tamanhos: 2048, 16384, 131072, 524288, 1048576 bp        │
│ • Formato: FASTA                                            │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. ANÁLISE FUNCIONAL COM ALPHAGENOME                        │
├─────────────────────────────────────────────────────────────┤
│ • Predição de expressão gênica                              │
│ • Acessibilidade cromatina                                  │
│ • Marcadores epigenéticos                                   │
│ • Fatores de transcrição                                    │
│ • 11 tipos de análise                                       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 5. CONSTRUÇÃO DO DATASET                                    │
├─────────────────────────────────────────────────────────────┤
│ Cada registro contém:                                       │
│ • Posição genômica (ponto central)                          │
│ • Sequência DNA (FASTA)                                     │
│ • Predições AlphaGenome                                     │
│ • Label (longevo=1, não-longevo=0)                          │
│ • Metadados                                                 │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 6. SPLITS: TREINO / VALIDAÇÃO / TESTE                       │
├─────────────────────────────────────────────────────────────┤
│ • 60% treino                                                │
│ • 20% validação                                             │
│ • 20% teste                                                 │
│ • Balanceamento de classes                                 │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 7. TREINAMENTO DE REDE NEURAL PROFUNDA                      │
├─────────────────────────────────────────────────────────────┤
│ • Arquitetura: a definir (CNN, Transformer, Híbrida)       │
│ • Framework: PyTorch                                        │
│ • Features: sequência DNA + predições AlphaGenome           │
│ • Loss: Binary Cross-Entropy                                │
│ • Métricas: Accuracy, Precision, Recall, F1, AUC-ROC       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│ 8. AVALIAÇÃO E INTERPRETAÇÃO                                │
├─────────────────────────────────────────────────────────────┤
│ • Teste em conjunto separado                                │
│ • Análise estatística de significância                      │
│ • Identificação de marcadores mais importantes             │
│ • Validação biológica                                       │
└─────────────────────────────────────────────────────────────┘
```

---

## 📁 Estrutura do Projeto

```
genomics/
├── neural_longevity_dataset.py     # Módulo principal
├── longevity_config.yaml           # Configuração
├── longevity_dataset/              # Dataset construído
│   ├── sequences/                  # Sequências FASTA
│   ├── alphagenome_cache/          # Cache de predições
│   ├── central_points.json         # Pontos centrais selecionados
│   ├── sequences_index.json        # Índice de sequências
│   ├── train.pkl                   # Dataset de treino (PyTorch)
│   ├── val.pkl                     # Dataset de validação
│   ├── test.pkl                    # Dataset de teste
│   └── checkpoint.json             # Checkpoint de processamento
├── neural_module.py                # Interface AlphaGenome
└── genomes_analyzer.py             # Pipeline genômico
```

---

## ⚙️ Configuração

### Arquivo: `longevity_config.yaml`

```yaml
project:
  name: "longevity_markers"
  output_dir: "longevity_dataset"

data_sources:
  longevous:
    source: "1000genomes_30x_grch38"
    sample_range: [0, 30]  # Primeiras 30 amostras
    n_samples: 30
    label: 1
  
  non_longevous:
    source: "1000genomes_30x_grch38"
    sample_range: [500, 530]  # 30 amostras a partir da posição 500
    n_samples: 30
    label: 0

variant_selection:
  initial_strategy: "first_longevous_sample"
  random_seed: 42
  n_central_points: 10
  filters:
    min_quality: 30
    min_depth: 10
    filter_pass_only: true

sequence_extraction:
  window_size: 2048  # Tamanho suportado pelo AlphaGenome
  use_alternate_allele: true

alphagenome:
  api_key: "YOUR_API_KEY"
  outputs:
    - "RNA_SEQ"
    - "CAGE"
    - "ATAC"
    - "CHIP_HISTONE"
    - "CHIP_TF"

dataset:
  splits:
    train: 0.6
    validation: 0.2
    test: 0.2
  balance_classes: true
  random_seed: 42
```

---

### Estratégias de Seleção de Pontos Centrais

- `first_longevous_sample`: usa apenas a primeira amostra longeva disponível como fonte para ordenar as variantes por QUAL e selecionar as `n_central_points` mais altas.
- `random_rotation_longevous_samples`: rotaciona todas as amostras longevas disponíveis. A cada iteração escolhe um longevo aleatório, seleciona uma variante elegível desse indivíduo **sem reposição** e incrementa o contador até atingir `n_central_points`. Caso alguma amostra não possua VCF ou variantes válidas, ela é descartada da rotação. O campo opcional `variant_selection.random_seed` garante reprodutibilidade dessa estratégia.

---

## 🚀 Como Usar

### 1. Preparação do Ambiente

```bash
conda activate genomics

# Instalar dependências adicionais
pip install torch torchvision torchaudio
pip install scikit-learn
```

### 2. Configurar API Key do AlphaGenome

Edite `longevity_config.yaml` e adicione sua API key:

```yaml
alphagenome:
  api_key: "YOUR_ALPHAGENOME_API_KEY"
```

### 3. Construir o Dataset

```bash
# Pipeline completo
python neural_longevity_dataset.py --config longevity_config.yaml

# Apenas etapas específicas
python neural_longevity_dataset.py \
  --config longevity_config.yaml \
  --steps download_samples select_central_points

# Modo simulação (dry-run)
python neural_longevity_dataset.py \
  --config longevity_config.yaml \
  --dry-run
```

### 4. Usar o Dataset em PyTorch

```python
import torch
from torch.utils.data import DataLoader
from neural_longevity_dataset import LongevityDataset

# Carregar datasets
train_dataset = LongevityDataset('longevity_dataset/train.pkl')
val_dataset = LongevityDataset('longevity_dataset/val.pkl')
test_dataset = LongevityDataset('longevity_dataset/test.pkl')

# Criar DataLoaders
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

# Iterar sobre dados
for features, labels in train_loader:
    sequence = features['sequence']      # DNA one-hot encoded
    position = features['position']      # Posição normalizada
    alphagenome = features['alphagenome'] # Features AlphaGenome
    
    # Seu modelo aqui
    # output = model(sequence, position, alphagenome)
    # loss = criterion(output, labels)
```

---

## 📊 Estrutura do Dataset

### Cada Amostra Contém:

1. **Sequência DNA** (one-hot encoded)
   - Shape: `(4, window_size)` para window_size=2048
   - Encoding: A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1]

2. **Posição Genômica** (normalizada)
   - Cromossomo e posição normalizados para [0, 1]

3. **Predições AlphaGenome** (features agregadas)
   - Para cada output (RNA_SEQ, CAGE, etc.):
     - mean, std, min, max, median
   - Total: ~35 features numéricas

4. **Label**
   - 1 = Longevo
   - 0 = Não-longevo

5. **Metadados**
   - `sample_id`: ID da amostra
   - `chromosome`: Cromossomo
   - `position`: Posição exata
   - `ref_allele`: Alelo de referência
   - `alt_allele`: Alelo alternativo
   - `variant_type`: Tipo de variante (SNV, INSERTION, DELETION)

---

## 🧠 Estratégia de Treinamento

### Fase 1: Exploração Inicial (ATUAL)

- **Dados**: Simulação com 1000 Genomes
- **Amostras**: 30 longevas + 30 não-longevas
- **Pontos**: 10 variantes da primeira pessoa longeva
- **Objetivo**: Validar pipeline e estrutura

### Fase 2: Seleção de Marcadores

Após primeira análise:
1. Treinar modelo inicial
2. Analisar importância de features
3. Identificar top-K pontos centrais mais relevantes
4. Reconstruir dataset com pontos refinados

### Fase 3: Dados Reais de Longevidade

Quando disponíveis:
1. Substituir dados simulados por dados reais
2. Aumentar número de amostras
3. Refinar seleção de pontos centrais
4. Treinar modelo final

---

## 🔬 Análise e Refinamento

### Identificação de Marcadores Importantes

```python
# TODO: Implementar análise de importância
# - Attention weights (se usar Transformers)
# - Gradient-based attribution
# - SHAP values
# - Permutation importance
```

### Reconstrução do Dataset

Após identificar os K pontos mais relevantes:

```bash
# Atualizar config com novos pontos centrais
# Reexecutar pipeline
python neural_longevity_dataset.py \
  --config longevity_config_refined.yaml \
  --steps select_central_points extract_sequences run_alphagenome build_dataset
```

---

## 📈 Arquiteturas de Rede Sugeridas

### Opção 1: CNN 1D

Adequada para capturar padrões locais na sequência DNA.

```python
class LongevityCNN(nn.Module):
    def __init__(self):
        super().__init__()
        # DNA sequence processing
        self.conv1 = nn.Conv1d(4, 64, kernel_size=7)
        self.conv2 = nn.Conv1d(64, 128, kernel_size=5)
        
        # AlphaGenome features processing
        self.fc_alpha = nn.Linear(35, 128)
        
        # Fusion and classification
        self.fc1 = nn.Linear(128 + 128, 256)
        self.fc2 = nn.Linear(256, 1)
    
    def forward(self, sequence, position, alphagenome):
        # Process sequence
        x_seq = F.relu(self.conv1(sequence))
        x_seq = F.max_pool1d(x_seq, 2)
        x_seq = F.relu(self.conv2(x_seq))
        x_seq = F.adaptive_max_pool1d(x_seq, 1).squeeze(-1)
        
        # Process AlphaGenome features
        x_alpha = F.relu(self.fc_alpha(alphagenome))
        
        # Fuse
        x = torch.cat([x_seq, x_alpha], dim=1)
        x = F.relu(self.fc1(x))
        x = torch.sigmoid(self.fc2(x))
        return x
```

### Opção 2: Transformer

Adequado para capturar dependências de longo alcance.

```python
class LongevityTransformer(nn.Module):
    # TODO: Implementar arquitetura Transformer
    pass
```

### Opção 3: Híbrida (CNN + Attention)

Combina vantagens de ambas abordagens.

---

## 📊 Métricas de Avaliação

### Métricas Principais

- **Accuracy**: Proporção de predições corretas
- **Precision**: TP / (TP + FP)
- **Recall (Sensitivity)**: TP / (TP + FN)
- **F1-Score**: Média harmônica de Precision e Recall
- **AUC-ROC**: Área sob curva ROC
- **Confusion Matrix**: Matriz de confusão

### Análise Estatística

- **Teste t**: Comparar performance vs baseline
- **Intervalo de confiança**: Bootstrap 95%
- **Cross-validation**: K-fold para robustez
- **Significância**: p-value < 0.05

---

## 🔄 Workflow de Refinamento

```
1. Treinar modelo inicial
        ↓
2. Avaliar no conjunto de validação
        ↓
3. Analisar importância de features
        ↓
4. Identificar top-K pontos centrais mais discriminativos
        ↓
5. Reconstruir dataset com pontos refinados
        ↓
6. Re-treinar modelo
        ↓
7. Avaliar no conjunto de teste (apenas no final)
        ↓
8. Validação biológica dos marcadores encontrados
```

---

## 🛠️ Próximos Passos

### Implementação Pendente

1. ✅ Estrutura básica do dataset builder
2. ✅ Configuração YAML
3. ✅ Classes de dados (GenomicVariant, SequenceRecord, etc.)
4. ✅ PyTorch Dataset class
5. ⏳ Completar extração de features do AlphaGenome
6. ⏳ Implementar build completo do dataset PyTorch
7. ⏳ Criar script de treinamento
8. ⏳ Implementar análise de importância
9. ⏳ Pipeline de refinamento automático
10. ⏳ Documentação de modelos

### Validação

1. Testar com dados simulados
2. Validar pipeline end-to-end
3. Benchmark de performance
4. Testes de robustez

### Dados Reais

1. Identificar fontes de genomas de pessoas longevas
2. Negociar acesso aos dados
3. Validação ética e LGPD
4. Processamento em larga escala

---

## 📚 Referências

### Bioinformática

- 1000 Genomes Project: https://www.internationalgenome.org/
- GRCh38 Reference: https://www.ncbi.nlm.nih.gov/grc/human
- VCF Format: https://samtools.github.io/hts-specs/VCFv4.3.pdf

### AlphaGenome

- Documentation: https://www.alphagenomedocs.com/
- Paper: (a ser publicado)

### Deep Learning

- PyTorch: https://pytorch.org/
- PyTorch Geometric: https://pytorch-geometric.readthedocs.io/
- Transformers: https://huggingface.co/docs/transformers/

### Longevidade

- Longevity Genes Database: https://genomics.senescence.info/
- Human Ageing Genomic Resources: https://genomics.senescence.info/

---

## 📝 Notas Importantes

### Limitações Atuais

1. **Dados Simulados**: Usando 1000 Genomes como proxy
2. **Pontos Centrais**: Seleção inicial arbitrária (top-10 por qualidade)
3. **AlphaGenome**: Custo de API (gratuito mas com limites)
4. **Dataset Small**: Apenas 60 amostras inicialmente

### Considerações Éticas

- Privacidade dos dados genômicos
- Consentimento informado
- LGPD e GDPR compliance
- Uso responsável de IA em saúde
- Validação clínica necessária antes de aplicação

### Performance

- Processamento intensivo (AlphaGenome API)
- Cache essencial para reutilização
- Paralelização quando possível
- Estimativa: ~5-10 min por amostra com AlphaGenome

---

**Status**: 🚧 EM DESENVOLVIMENTO

**Última atualização**: Outubro 2025

**Contato**: Projeto IA Neuro-Simbólica para Longevidade

