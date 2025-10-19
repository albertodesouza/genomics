# 💡 Guia de Uso - Neural Module

## 🎯 Visão Geral

Este guia mostra como usar o Neural Module para análise de DNA com AlphaGenome.

**💡 Primeiro Uso?** Veja como baixar sequências genômicas reais: **[Download de Sequências](DOWNLOAD_SEQUENCES.md)**

---

## 🚀 Uso Básico

### Comando Mínimo
```bash
python neural_module.py -i INPUT.fasta -k API_KEY -o OUTPUT_DIR/
```

### Exemplo Real
```bash
python neural_module.py \
    -i example_sequence.fasta \
    -k AIzaSyAUFo... \
    -o my_results/
```

---

## 📊 Tipos de Análise

### 1. Análise de Sequência (Padrão)

Analisa características funcionais de uma sequência de DNA.

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/
```

**Saídas**:
- Gráficos de RNA-seq, CAGE, ATAC, CHIP_HISTONE, CHIP_TF
- Heatmaps (múltiplas tracks)
- Dashboard resumo
- Comparação entre outputs
- analysis_report.json

---

### 2. Análise de Variante

Prediz o efeito de uma mutação específica.

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --variant POSIÇÃO REF ALT
```

**Exemplo - Anemia Falciforme**:
```bash
python neural_module.py \
    -i example_sequence.fasta \
    -k API_KEY \
    -o sickle_cell/ \
    --variant 1024 A T
```

**Saídas**:
- 3 painéis: Sobreposição REF vs ALT, Diferença, Zoom
- Marcação da posição da variante
- Gráficos comparativos

---

### 3. Análise com Outputs Específicos

Escolha quais tipos de análise executar:

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ ATAC
```

**Outputs disponíveis**:
- `RNA_SEQ` - Expressão gênica via RNA-seq
- `CAGE` - Cap Analysis of Gene Expression
- `ATAC` - Acessibilidade de cromatina
- `CHIP_HISTONE` - Marcadores de histonas
- `CHIP_TF` - Fatores de transcrição (CTCF, etc.)
- `DNASE` - DNase-seq
- `PROCAP` - PRO-cap
- `CONTACT_MAPS` - Mapas de contato 3D
- `SPLICE_JUNCTIONS` - Junções de splicing
- `SPLICE_SITES` - Sítios de splicing
- `SPLICE_SITE_USAGE` - Uso de sítios

Ver todos: [OUTPUTS_DISPONIVEIS.md](OUTPUTS_DISPONIVEIS.md)

---

### 4. Alta Resolução para Publicação

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o publication/ \
    --dpi 600 \
    --formats png pdf svg
```

**Configurações**:
- `--dpi 600` - Alta resolução (padrão: 300)
- `--formats png pdf svg` - Múltiplos formatos

---

### 5. Análise com Metadados de Ontologia

Salva informações detalhadas sobre tecidos/células analisados (HABILITADO POR PADRÃO):

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/
    # Metadados são salvos automaticamente!
```

**Saídas adicionais**:
- `seq_id_RNA_SEQ_metadata.csv` - Metadados em CSV
- `seq_id_RNA_SEQ_metadata.json` - Metadados em JSON
- Tabela no terminal com informações dos tecidos

Para **desabilitar** metadados:
```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --no-metadata
```

### 6. Análise Rápida (Sem Visualizações)

Para análises rápidas ou quando só precisa dos dados:

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --no-plots
```

**Saída**: Apenas `analysis_report.json` e metadados

---

## ⚙️ Opções Avançadas

### Especificar Contexto Genômico

```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --chromosome chr11 \
    --start 5246696
```

Útil para:
- Referências específicas
- Análise de regiões conhecidas
- Correlação com dados externos

---

### Controlar Visualizações

As visualizações avançadas são **padrão**. Incluem:
- ✅ Enhanced tracks (múltiplos subplots)
- ✅ Heatmaps
- ✅ Dashboard resumo
- ✅ Comparação multi-output

Visualizações são sempre habilitadas por padrão.

---

## 📁 Formato de Entrada (FASTA)

### Requisitos

1. **Formato FASTA padrão**:
```fasta
>sequence_id description
ATCGATCGATCG...
```

2. **Tamanhos suportados** (CRÍTICO!):
   - ✅ 2,048 bp (2 KB)
   - ✅ 16,384 bp (16 KB)
   - ✅ 131,072 bp (128 KB)
   - ✅ 524,288 bp (512 KB)
   - ✅ 1,048,576 bp (1 MB)

**Qualquer outro tamanho resultará em erro!**

Ver detalhes: [TAMANHOS_SUPORTADOS.md](TAMANHOS_SUPORTADOS.md)

---

## 📊 Entendendo as Saídas

### Arquivos Gerados

```
results/
├── seq_id_RNA_SEQ.png              # Plot básico RNA-seq
├── seq_id_RNA_SEQ_enhanced.png     # Enhanced (subplots + metadados)
├── seq_id_RNA_SEQ_heatmap.png      # Heatmap de tracks
├── seq_id_CAGE.png
├── seq_id_CAGE_enhanced.png
├── seq_id_ATAC.png
├── seq_id_ATAC_enhanced.png
├── seq_id_comparison.png           # Comparação multi-output
├── seq_id_dashboard.png            # Dashboard resumo
└── analysis_report.json            # Relatório JSON
```

### Interpretando Gráficos

Ver guia detalhado: **[RESULTS_NEURAL.md](RESULTS_NEURAL.md)**

---

## 💻 Uso Programático

### Como Biblioteca Python

```python
from neural_module import AlphaGenomeAnalyzer

# Inicializar
analyzer = AlphaGenomeAnalyzer(api_key="YOUR_KEY")
analyzer.initialize()

# Analisar sequência
resultado = analyzer.predict_sequence(
    sequence="ATCG" * 512,  # 2048 bp
    seq_id="minha_seq",
    requested_outputs=["RNA_SEQ", "ATAC"]
)

# Analisar variante
resultado_var = analyzer.predict_variant(
    sequence="ATCG" * 512,
    seq_id="variante",
    variant_position=1000,
    ref_base="A",
    alt_base="T"
)
```

Ver exemplos completos: [neural_example.py](neural_example.py)

---

## 🔗 Integração com Pipeline

Para integrar o Neural Module com o `genomes_analyzer.py` e analisar variantes automaticamente:

### Análise Integrada Completa

```bash
python neural_integration.py \
    --integrated \
    --vcf vcf/sample.vcf.gz \
    --ref refs/GRCh38.fa \
    --api-key YOUR_API_KEY \
    --output integrated_results/
```

### Extrair Sequências de VCF

```bash
python neural_integration.py \
    --extract-vcf \
    --vcf variants.vcf \
    --ref genome.fa \
    --output extracted.fasta
```

### Extrair Genes Específicos

```bash
python neural_integration.py \
    --extract-genes \
    --genes BRCA1 TP53 HBB \
    --gtf refs/gencode.gtf.gz \
    --ref refs/GRCh38.fa \
    --output genes.fasta \
    --flank 10000
```

📖 **Guia Completo de Integração**: [NEURAL_INTEGRATION.md](NEURAL_INTEGRATION.md)

O guia de integração inclui:
- 4 modos de operação detalhados
- Casos de uso práticos
- Fluxo completo pipeline → neural
- Interpretação de resultados correlacionados

---

## 🧬 Exemplo Prático: Anemia Falciforme

### Contexto Biológico

A anemia falciforme é causada por uma mutação no gene **HBB** (beta-globina):
- **Posição**: chr11:5227002 (hg38)
- **Mutação**: A→T (GAG→GTG)
- **Efeito**: Glutamato → Valina (posição 6)
- **Consequência**: Hemoglobina defeituosa (HbS)

### Análise com Neural Module

```bash
# 1. Analisar região do gene HBB
python neural_module.py \
    -i example_sequence.fasta \
    -k API_KEY \
    -o hbb_analysis/

# 2. Analisar efeito da mutação
python neural_module.py \
    -i example_sequence.fasta \
    -k API_KEY \
    -o sickle_cell_variant/ \
    --variant 1024 A T

# 3. Comparar múltiplos outputs
python neural_module.py \
    -i example_sequence.fasta \
    -k API_KEY \
    -o comprehensive/ \
    --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE \
    --variant 1024 A T
```

### Interpretação

1. **RNA_SEQ**: Mostra impacto na expressão do gene
2. **ATAC**: Indica mudanças na acessibilidade
3. **CHIP_HISTONE**: Revela alterações epigenéticas
4. **Comparação REF vs ALT**: Visualiza diferença causada pela mutação

---

## ⚡ Dicas de Performance

### 1. Começar Pequeno
```bash
# Teste primeiro com RNA_SEQ
python neural_module.py -i seq.fasta -k API_KEY -o test/ --outputs RNA_SEQ
```

### 2. Usar `--no-plots` para Dados Apenas
```bash
python neural_module.py -i seq.fasta -k API_KEY -o data_only/ --no-plots
```

### 3. Escolher Tamanho Adequado

- Variantes: 2KB ou 16KB
- Genes: 16KB ou 128KB
- Regiões grandes: 512KB ou 1MB

---

## 🐛 Troubleshooting

### Erro: "Sequence length X not supported"
**Solução**: Use um dos 5 tamanhos suportados (2KB, 16KB, 128KB, 512KB, 1MB)

### Erro: "Output 'X' not available"
**Solução**: Verifique outputs disponíveis:
```bash
python check_alphagenome_outputs.py
```

### Análise muito lenta
**Solução**: 
- Use sequências menores
- Reduza número de outputs
- Use `--no-plots` para análise rápida

---

## 📚 Recursos Adicionais

- **[README Completo](NEURAL_MODULE_README.md)** - Documentação técnica
- **[Outputs Disponíveis](OUTPUTS_DISPONIVEIS.md)** - Lista completa
- **[Visualizações Avançadas](VISUALIZACOES_AVANCADAS.md)** - Guia visual
- **[Tamanhos Suportados](TAMANHOS_SUPORTADOS.md)** - Detalhes técnicos

---

## ✅ Checklist de Uso

Antes de cada análise:

- [ ] Arquivo FASTA tem tamanho suportado?
- [ ] API key está correta?
- [ ] Outputs escolhidos são válidos?
- [ ] Diretório de saída tem espaço?
- [ ] Conexão com internet funcionando?

---

**Próximo passo**: [Interpretar Resultados →](RESULTS_NEURAL.md)

*Última atualização: Outubro 2025*

