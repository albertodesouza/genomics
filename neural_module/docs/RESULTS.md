# 📊 Interpretando Resultados - Neural Module

## 🎯 Visão Geral

Este guia explica como interpretar as visualizações e dados gerados pelo Neural Module.

---

## 📁 Estrutura de Saída

Após executar uma análise, você terá:

```
results/
├── seq_id_RNA_SEQ.png              # Plot básico
├── seq_id_RNA_SEQ_enhanced.png     # Com subplots e metadados
├── seq_id_RNA_SEQ_heatmap.png      # Heatmap de múltiplas tracks
├── seq_id_CAGE.png
├── seq_id_CAGE_enhanced.png
├── seq_id_ATAC.png
├── seq_id_ATAC_enhanced.png
├── seq_id_CHIP_HISTONE.png
├── seq_id_CHIP_HISTONE_enhanced.png
├── seq_id_CHIP_HISTONE_heatmap.png
├── seq_id_CHIP_TF.png
├── seq_id_CHIP_TF_enhanced.png
├── seq_id_comparison.png           # 🆕 Comparação entre outputs
├── seq_id_dashboard.png            # 🆕 Dashboard resumo
└── analysis_report.json            # Dados estruturados
```

---

## 🎨 Tipos de Visualização

### 1️⃣ Plot Básico

**Arquivo**: `*_OUTPUT.png`

**O que mostra**:
- Sinal ao longo da sequência
- Eixo X: Posição (bp)
- Eixo Y: Intensidade do sinal

**Como interpretar**:
- **Picos** = Regiões ativas
- **Vales** = Regiões inativas
- **Padrões** = Estruturas funcionais

---

### 2️⃣ Enhanced (Melhorado)

**Arquivo**: `*_OUTPUT_enhanced.png`

**O que mostra**:
- Múltiplos subplots (um por track/tecido)
- Metadados visíveis:
  - Tecido (ex: "Brain", "Liver")
  - Strand (+ ou -)
  - Assay type

**Como interpreter**:
```
┌─────────────────────────────────────┐
│ Tecido: Brain | Strand: + | RNA-seq │
│ ▁▂▃▄▅▆▇█  ← Expressão alta no cérebro
├─────────────────────────────────────┤
│ Tecido: Liver | Strand: + | RNA-seq │
│ ▁▁▂▂▃▃    ← Expressão baixa no fígado
└─────────────────────────────────────┘
```

**Vantagens**:
- ✅ Compara diferentes tecidos
- ✅ Identifica especificidade tecidual
- ✅ Mostra strand-specificity

---

### 3️⃣ Heatmap

**Arquivo**: `*_OUTPUT_heatmap.png`

**O que mostra**:
- Matriz de múltiplas tracks
- Cor = Intensidade do sinal
- Linhas = Tecidos/condições
- Colunas = Posições na sequência

**Como interpretar**:
```
        Posição (bp) →
Brain   █▓▒░░▒▓█  ← Padrão de expressão
Liver   ▓▒░░▒▓█▓
Lung    ░░▒▓█▓▒░
```

- **Vermelho/Amarelo** = Sinal alto
- **Azul/Roxo** = Sinal baixo
- **Padrões verticais** = Regiões funcionais conservadas
- **Padrões horizontais** = Especificidade tecidual

---

### 4️⃣ Comparação Multi-Output

**Arquivo**: `*_comparison.png`

**O que mostra**:
- Todos outputs em um gráfico
- Subplots empilhados
- Estatísticas (μ, σ) por output

**Como interpretar**:
```
RNA_SEQ  ▄▅▆▇█▇▆▅▄  μ=0.123 σ=0.045
ATAC     ▁▂▃▄▅▄▃▂▁  μ=0.234 σ=0.089
CHIP_H   ▂▃▄▅▆▅▄▃▂  μ=0.156 σ=0.067
```

**Insights**:
- **Correlação** entre outputs?
- **Regiões ativas** em múltiplas dimensões?
- **Hotspots regulatórios**?

---

### 5️⃣ Dashboard

**Arquivo**: `*_dashboard.png`

**O que mostra**:
- 4 painéis:
  1. Comparação de médias (bar chart)
  2. Número de tracks por output
  3. Range de valores
  4. Tabela de estatísticas

**Como interpretar**:

```
┌──────────────────┬──────────────────┐
│ Comparação       │  N Tracks        │
│ ▆ RNA ▄ ATAC    │  RNA    █████ 5  │
├──────────────────┴──────────────────┤
│ Estatísticas                        │
│ Out │ μ    │ σ    │ min │ max     │
│ RNA │0.123 │0.045 │0.01 │0.87     │
└─────────────────────────────────────┘
```

**Uso**:
- ✅ Resumo rápido da análise
- ✅ Identificar outputs mais informativos
- ✅ Quality control

---

## 🧬 Análise de Variantes

### Arquivo: `*_variant_enhanced.png`

**3 Painéis**:

#### 1. Sobreposição REF vs ALT
```
REF: ___/‾‾‾\___  (cinza)
ALT: ___/‾‾‾\___  (vermelho)
       ↑ variante
```

**Interpretar**:
- **Coincidência** = Sem efeito?
- **Divergência** = Efeito funcional
- **Amplitude** = Magnitude do efeito

#### 2. Diferença (ALT - REF)
```
Δ:  ___/‾‾\___ 
      ↑ aumento (verde)
      ↓ diminuição (vermelho)
```

**Interpretar**:
- **Verde** = Variante **aumenta** sinal
- **Vermelho** = Variante **diminui** sinal
- **Magnitude** = Intensidade do efeito

#### 3. Zoom (±200bp)
```
Zoom da região ao redor da variante
REF: •—•—•—|—•—•—•
ALT: •—•—•—|—•—•—•
           ↑ variante
```

**Interpretar**:
- Efeito **local** vs **distal**
- **Propagação** do efeito
- **Padrões** finos

---

## 📈 Interpretação por Output Type

### RNA_SEQ (Expressão Gênica)

**O que indica**:
- Nível de transcrição do gene
- Atividade transcricional

**Picos** = Genes ativos  
**Vales** = Genes silenciados

**Variantes**:
- ↑ Sinal = Aumento de expressão
- ↓ Sinal = Diminuição de expressão

---

### CAGE (Cap Analysis)

**O que indica**:
- Início de transcrição (TSS)
- Promotores ativos

**Picos afiados** = Promotores  
**Múltiplos picos** = TSSs alternativos

**Variantes**:
- Mudança em pico = Mudança em promotor
- Novo pico = Novo TSS

---

### ATAC (Acessibilidade)

**O que indica**:
- Cromatina aberta/acessível
- Regiões regulatórias ativas

**Picos** = Enhancers, promotores  
**Padrão periódico** = Nucleossomos

**Variantes**:
- ↑ Acessibilidade = Ativação
- ↓ Acessibilidade = Repressão

---

### CHIP_HISTONE (Marcadores)

**O que indica**:
- Estado epigenético
- Múltiplas marcas simultâneas

**Tipos de marcadores**:
- H3K27AC = Enhancers ativos
- H3K4ME3 = Promotores ativos
- H3K27ME3 = Repressão (Polycomb)
- H3K36ME3 = Corpo gênico ativo
- H3K9ME3 = Heterocromatina

**Variantes**:
- Mudança em marca = Mudança epigenética

---

### CHIP_TF (Fatores de Transcrição)

**O que indica**:
- Sítios de ligação de TFs
- Regulação transcricional

**CTCF**:
- Picos = Insuladores
- Organizadores de loops de cromatina

**Variantes**:
- Ganho/perda de pico = Ganho/perda de binding site

---

## 📊 Relatório JSON

**Arquivo**: `analysis_report.json`

```json
{
  "timestamp": "2025-10-18T15:00:00",
  "total_sequences": 1,
  "successful_analyses": 1,
  "sequences": [
    {
      "id": "HBB_gene",
      "length": 2048,
      "status": "success",
      "outputs": ["RNA_SEQ", "CAGE", "ATAC", "CHIP_HISTONE", "CHIP_TF"]
    }
  ]
}
```

**Uso**:
- Processamento automatizado
- Integração com pipelines
- Quality control

---

## 🎓 Caso de Estudo: Anemia Falciforme

### Contexto

Mutação no gene HBB: A→T (posição 1024)
- **Normal**: GAG (Glutamato)
- **Mutado**: GTG (Valina)

### Análise

```bash
python neural_module.py \
    -i example_sequence.fasta \
    -k API_KEY \
    -o sickle_cell/ \
    --variant 1024 A T
```

### Interpretação Esperada

#### 1. RNA_SEQ
- **Diferença mínima**: Mutação não afeta muito a transcrição
- Explicação: Mudança é no coding, não regulatória

#### 2. ATAC
- **Sem mudança significativa**: Acessibilidade mantida
- Explicação: Mutação não cria/remove binding sites

#### 3. CHIP_HISTONE
- **Padrão similar**: Estado epigenético conservado
- Explicação: Mutação pontual tem efeito local

#### 4. Zoom na variante
- **Pequena alteração local**: Efeito estrutural da proteína
- Explicação: Impacto é na sequência de aminoácidos, não na regulação

### Conclusão

Anemia falciforme é causada por mudança **estrutural** da proteína, não mudança **regulatória**. Por isso, Neural Module mostra:
- ✅ Gene continua expresso
- ✅ Regulação mantida
- ⚠️ Mas proteína é defeituosa (não capturado pelo modelo)

---

## 💡 Dicas de Interpretação

### 1. Compare Múltiplos Outputs
```
RNA_SEQ alto + ATAC alto = Gene ativo e acessível ✓
RNA_SEQ baixo + ATAC alto = Regulação pós-transcricional?
```

### 2. Contexto Biológico
- Sabe qual tecido/célula interessa? Foque nessas tracks
- Conhece a função do gene? Espera certos padrões

### 3. Variantes
- **Grande efeito** = Provável patogenicidade
- **Sem efeito** = Provável benigna ou non-coding
- **Efeito tecido-específico** = Fenótipo restrito

### 4. Use Dashboard
- Visão geral rápida
- Identifique outliers
- Quality control

---

## 🔍 Análise Avançada

### Correlação Entre Outputs

Se RNA_SEQ e H3K4ME3 estão altos juntos:
→ Promotor ativo com transcrição

Se ATAC alto mas RNA_SEQ baixo:
→ Enhancer poised ou regulação pós-transcricional

Se H3K27ME3 alto:
→ Repressão Polycomb, gene silenciado

### Padrões Espaciais

**Clusters de picos** = Enhancer clusters, super-enhancers  
**Picos isolados** = Elementos regulatórios individuais  
**Padrão periódico** = Estrutura nucleossomal

---

## 📚 Recursos Adicionais

- **[Documentação AlphaGenome](https://www.alphagenomedocs.com/)** - Detalhes técnicos
- **[Paper](https://doi.org/10.1101/2025.06.25.661532)** - Metodologia
- **[OUTPUTS_DISPONIVEIS.md](OUTPUTS_DISPONIVEIS.md)** - Lista completa

---

## ✅ Checklist de Interpretação

- [ ] Vi todos os plots gerados?
- [ ] Verifiquei dashboard para QC?
- [ ] Comparei múltiplos outputs?
- [ ] Considerei contexto biológico?
- [ ] Para variantes: efeito faz sentido biologicamente?
- [ ] Consultei literature sobre o gene/região?

---

**Dúvidas?** Consulte:
- [Guia de Uso](USAGE_NEURAL.md)
- [README Completo](NEURAL_MODULE_README.md)
- [Documentação do AlphaGenome](https://www.alphagenomedocs.com/)

*Última atualização: Outubro 2025*

