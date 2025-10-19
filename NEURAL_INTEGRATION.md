# 🔗 Integração Neural Module + Genomes Analyzer

## 📖 Visão Geral

O **Neural Integration** (`neural_integration.py`) é uma ferramenta de ponte que conecta o pipeline tradicional de análise genômica (`genomes_analyzer.py`) com a análise neural baseada em IA (`neural_module.py` + AlphaGenome).

Ele automatiza o fluxo de:
**VCF/BED/Genes → Extração de Sequências → Análise Neural → Correlação de Resultados**

---

## 🎯 O que o Neural Integration faz?

### 1. **Extração Inteligente de Sequências**
- ✅ Extrai regiões do genoma a partir de VCF (variantes)
- ✅ Extrai regiões de arquivos BED (regiões de interesse)
- ✅ Extrai genes específicos de GTF (com regiões flanqueadoras)
- ✅ Converte tudo para FASTA pronto para AlphaGenome

### 2. **Análise Neural Automatizada**
- ✅ Executa `neural_module.py` automaticamente
- ✅ Configura parâmetros apropriados
- ✅ Gerencia outputs do AlphaGenome

### 3. **Correlação de Resultados**
- ✅ Correlaciona variantes com predições neurais
- ✅ Gera relatórios de integração
- ✅ Cria visualizações combinadas

---

## 🚀 Modos de Operação

O `neural_integration.py` opera em 4 modos diferentes:

### Modo 1: **Análise Integrada Completa** 🌟

Executa todo o fluxo automaticamente: extração + análise + correlação.

```bash
python neural_integration.py \
  --integrated \
  --vcf vcf/NA12878.vcf.gz \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output integrated_results/
```

**Quando usar**: Após executar `genomes_analyzer.py` e ter um VCF de variantes.

**O que acontece**:
1. Extrai regiões ao redor de cada variante (±5kb)
2. Converte para FASTA
3. Executa análise neural com AlphaGenome
4. Correlaciona variantes com predições
5. Gera relatórios integrados

---

### Modo 2: **Extrair Sequências de VCF**

Apenas extrai sequências sem executar análise neural.

```bash
python neural_integration.py \
  --extract-vcf \
  --vcf vcf/NA12878.vcf.gz \
  --ref refs/GRCh38.d1.vd1.fa \
  --output variants_sequences.fasta
```

**Quando usar**: Quando você quer apenas preparar sequências para análise posterior.

**Saída**: Arquivo FASTA com regiões de ±5kb ao redor de cada variante.

---

### Modo 3: **Extrair Sequências de BED**

Extrai sequências de regiões especificadas em arquivo BED.

```bash
python neural_integration.py \
  --extract-bed \
  --bed regions_of_interest.bed \
  --ref refs/GRCh38.d1.vd1.fa \
  --output regions_sequences.fasta
```

**Quando usar**: Quando você tem regiões de interesse específicas (enhancers, promotores, etc.).

**Formato BED esperado**:
```
chr1    1000000    1002048    region_1
chr2    5000000    5002048    region_2
```

---

### Modo 4: **Extrair Genes Específicos**

Extrai sequências de genes por nome, com regiões flanqueadoras.

```bash
python neural_integration.py \
  --extract-genes \
  --genes BRCA1 TP53 HBB CFTR \
  --gtf refs/gencode.v38.annotation.gtf.gz \
  --ref refs/GRCh38.d1.vd1.fa \
  --output genes_sequences.fasta \
  --flank 10000
```

**Quando usar**: Para analisar genes específicos de interesse com contexto regulatório.

**Parâmetros**:
- `--flank`: Bases flanqueadoras (padrão: 10kb antes e depois do gene)

---

## 📊 Casos de Uso Práticos

### Caso 1: Analisar Variantes de Alto Impacto

Após identificar variantes de alto impacto no VCF:

```bash
# Passo 1: Filtrar variantes de alto impacto (exemplo)
bcftools view -i 'INFO/ANN~"HIGH"' vcf/sample.vcf.gz > high_impact.vcf

# Passo 2: Análise integrada
python neural_integration.py \
  --integrated \
  --vcf high_impact.vcf \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output neural_high_impact/ \
  --outputs RNA_SEQ ATAC CHIP_HISTONE
```

**Resultado**: Predições neurais para regiões com variantes de alto impacto.

---

### Caso 2: Analisar Genes Candidatos

Você identificou genes candidatos de doença:

```bash
python neural_integration.py \
  --extract-genes \
  --genes BRCA1 BRCA2 TP53 PTEN \
  --gtf refs/gencode.v38.annotation.gtf.gz \
  --ref refs/GRCh38.d1.vd1.fa \
  --output candidate_genes.fasta \
  --flank 20000

# Depois analisar com neural_module
python neural_module.py \
  -i candidate_genes.fasta \
  -k YOUR_API_KEY \
  -o candidate_genes_neural/
```

**Resultado**: Análise funcional completa dos genes candidatos.

---

### Caso 3: Regiões Regulatórias Não-Codificantes

Você tem regiões regulatórias de interesse em BED:

```bash
# enhancers.bed contém regiões de enhancers
python neural_integration.py \
  --extract-bed \
  --bed enhancers.bed \
  --ref refs/GRCh38.d1.vd1.fa \
  --output enhancers.fasta

# Analisar acessibilidade e marcadores epigenéticos
python neural_module.py \
  -i enhancers.fasta \
  -k YOUR_API_KEY \
  -o enhancers_analysis/ \
  --outputs ATAC DNASE CHIP_HISTONE CHIP_TF
```

**Resultado**: Predições de atividade regulatória.

---

### Caso 4: Análise de Trio (Variantes De Novo)

Após identificar variantes de novo no trio:

```bash
# Passo 1: Obter variantes de novo do pipeline
# (genomes_analyzer.py já gera trio/denovo_candidates.vcf)

# Passo 2: Análise neural das variantes de novo
python neural_integration.py \
  --integrated \
  --vcf trio/denovo_candidates.vcf \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output neural_denovo/ \
  --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF
```

**Resultado**: Impacto funcional predito das variantes de novo.

---

## 🔄 Fluxo de Trabalho Integrado Completo

### Pipeline Completo: DNA → Variantes → Predições Neurais

```bash
# ════════════════════════════════════════════════════════════════
# PASSO 1: Análise Genômica Tradicional
# ════════════════════════════════════════════════════════════════

conda activate genomics
python genomes_analyzer.py --config config_human_30x.yaml

# Saídas:
# - vcf/NA12878.vcf.gz (variantes)
# - trio/denovo_candidates.vcf (de novo)
# - bam/*.bam (alinhamentos)

# ════════════════════════════════════════════════════════════════
# PASSO 2: Filtrar Variantes de Interesse (Opcional)
# ════════════════════════════════════════════════════════════════

# Exemplo: Variantes exônicas de alto impacto
bcftools view -i 'INFO/ANN~"HIGH|MODERATE" && INFO/ANN~"exonic"' \
  vcf/NA12878.vcf.gz > variants_of_interest.vcf

# ════════════════════════════════════════════════════════════════
# PASSO 3: Análise Neural Integrada
# ════════════════════════════════════════════════════════════════

python neural_integration.py \
  --integrated \
  --vcf variants_of_interest.vcf \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_ALPHAGENOME_KEY \
  --output integrated_analysis/ \
  --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF

# Saídas:
# - integrated_analysis/neural_results/ (predições AlphaGenome)
# - integrated_analysis/correlation_report.json (correlação)

# ════════════════════════════════════════════════════════════════
# PASSO 4: Interpretar Resultados
# ════════════════════════════════════════════════════════════════

# Ver relatório de correlação
cat integrated_analysis/correlation_report.json | jq .

# Ver visualizações neurais
ls integrated_analysis/neural_results/*.png

# Ver metadados de ontologia
cat integrated_analysis/neural_results/*_metadata.csv
```

---

## 📁 Estrutura de Saídas

### Modo Integrado (`--integrated`)

```
integrated_analysis/
├── variants_sequences.fasta          # Sequências extraídas do VCF
├── neural_results/                   # Resultados do neural_module
│   ├── variant_1_*_RNA_SEQ.png
│   ├── variant_1_*_RNA_SEQ_enhanced.png
│   ├── variant_1_*_RNA_SEQ_heatmap.png
│   ├── variant_1_*_RNA_SEQ_metadata.csv
│   ├── variant_1_*_RNA_SEQ_metadata.json
│   ├── ... (outros outputs)
│   ├── variant_1_*_comparison.png
│   ├── variant_1_*_dashboard.png
│   └── analysis_report.json
└── correlation_report.json           # Correlação variantes × predições
```

### Modo Extração (`--extract-*`)

```
sequences.fasta                       # Sequências extraídas prontas para uso
```

---

## 💡 Exemplos Avançados

### Exemplo 1: Análise por Cromossomo

```bash
# Analisar apenas variantes do chr11 (gene HBB)
python neural_integration.py \
  --integrated \
  --vcf <(bcftools view -r chr11 vcf/sample.vcf.gz) \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output chr11_neural/
```

### Exemplo 2: Genes de uma Via Metabólica

```bash
# Extrair todos os genes de uma via (exemplo: reparo de DNA)
python neural_integration.py \
  --extract-genes \
  --genes BRCA1 BRCA2 ATM CHEK2 TP53 PALB2 RAD51 \
  --gtf refs/gencode.v38.annotation.gtf.gz \
  --ref refs/GRCh38.d1.vd1.fa \
  --output dna_repair_genes.fasta \
  --flank 15000
```

### Exemplo 3: Priorização de Variantes

```bash
#!/bin/bash
# Script para priorizar variantes com análise neural

# 1. Variantes raras de alto impacto
bcftools view -i 'INFO/AF<0.01 && INFO/ANN~"HIGH"' \
  vcf/sample.vcf.gz > rare_high_impact.vcf

# 2. Análise neural
python neural_integration.py \
  --integrated \
  --vcf rare_high_impact.vcf \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key $ALPHAGENOME_KEY \
  --output prioritized_variants/ \
  --outputs RNA_SEQ ATAC CHIP_HISTONE

# 3. Ver resultados
python -c "
import json
with open('prioritized_variants/correlation_report.json') as f:
    data = json.load(f)
    print(f'Variantes analisadas: {data[\"summary\"][\"total_sequences\"]}')
    print(f'Predições bem-sucedidas: {data[\"summary\"][\"successful_predictions\"]}')
"
```

---

## ⚙️ Configuração Avançada

### Personalizar Tamanho de Janela

Por padrão, extrai ±5kb ao redor de variantes. Para mudar, edite `neural_integration.py`:

```python
# Linha ~76
start = max(1, pos - 10000)  # Era 5000, agora 10kb
end = pos + 10000
```

### Adicionar Filtros de Qualidade

Filtrar variantes antes de análise neural:

```bash
# Apenas variantes PASS com DP ≥ 20 e GQ ≥ 30
bcftools view -f PASS -i 'FORMAT/DP>=20 && FORMAT/GQ>=30' \
  vcf/sample.vcf.gz | \
python neural_integration.py \
  --integrated \
  --vcf /dev/stdin \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output high_quality_neural/
```

---

## 🔍 Interpretando Resultados

### Relatório de Correlação

O arquivo `correlation_report.json` contém:

```json
{
  "vcf_source": "vcf/sample.vcf.gz",
  "neural_results": {
    "timestamp": "2025-10-18T...",
    "total_sequences": 42,
    "successful_analyses": 40,
    "sequences": [
      {
        "id": "variant_1_chr11_5227002",
        "length": 10000,
        "status": "success",
        "outputs": ["RNA_SEQ", "CAGE", "ATAC", "CHIP_HISTONE", "CHIP_TF"]
      },
      ...
    ]
  },
  "summary": {
    "total_sequences": 42,
    "successful_predictions": 40
  }
}
```

### Análise de Resultados

1. **Visualizações Individuais**: Veja `neural_results/*_enhanced.png` para cada variante
2. **Heatmaps**: Compare múltiplas tracks em `*_heatmap.png`
3. **Dashboard**: Resumo estatístico em `*_dashboard.png`
4. **Metadados**: Informações de tecidos/células em `*_metadata.csv`

---

## ❓ FAQ

### P: O neural_integration requer o genomes_analyzer instalado?

**R**: Não! O `neural_integration.py` é independente. Ele apenas requer:
- `bcftools` (para VCF)
- `bedtools` (para BED)
- `samtools` (para extração de sequências)
- `neural_module.py` (para análise neural)

Todos já estão no ambiente `genomics`.

### P: Posso usar com VCFs de outros pipelines?

**R**: Sim! Funciona com qualquer VCF padrão, não apenas os gerados por `genomes_analyzer.py`.

### P: Quanto custa usar com AlphaGenome?

**R**: AlphaGenome é gratuito para uso não comercial. Veja https://www.alphagenomedocs.com/

### P: Posso analisar apenas algumas variantes específicas?

**R**: Sim! Use `bcftools view` para filtrar o VCF primeiro:

```bash
# Apenas variantes em posições específicas
bcftools view -t chr11:5227002 vcf/sample.vcf.gz | \
python neural_integration.py --integrated --vcf /dev/stdin ...
```

### P: Como adicionar análise de variantes (REF vs ALT)?

**R**: Use `neural_module.py` diretamente após extração:

```bash
# 1. Extrair sequências
python neural_integration.py \
  --extract-vcf \
  --vcf variants.vcf \
  --ref genome.fa \
  --output sequences.fasta

# 2. Analisar cada variante com --variant
# (requer script adicional para parse do VCF)
```

---

## 🔗 Recursos Relacionados

- **[Neural Module Principal](NEURAL_MODULE.md)** - Documentação completa
- **[Guia de Uso](USAGE_NEURAL.md)** - Como usar neural_module.py
- **[Download de Sequências](DOWNLOAD_SEQUENCES.md)** - Baixar genomas reais
- **[Interpretação de Resultados](RESULTS_NEURAL.md)** - Entender predições

---

## 🚀 Próximos Passos

Após dominar a integração básica:

1. **Automatizar**: Crie scripts para análise em lote
2. **Priorizar**: Combine scores de variantes com predições neurais
3. **Validar**: Compare predições com dados experimentais (se disponível)
4. **Publicar**: Inclua análises neurais em seus relatórios

---

**Criado**: Outubro 2025  
**Parte do**: Neural Module Documentation

