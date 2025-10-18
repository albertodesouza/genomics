# 📊 Outputs Disponíveis no AlphaGenome

## Lista Completa de Outputs

O AlphaGenome oferece **11 tipos de outputs** diferentes para análise de DNA:

### 🧬 Expressão Gênica

#### 1. `RNA_SEQ`
- **Descrição**: Predição de níveis de expressão gênica via RNA-seq
- **Uso**: Identificar genes ativos e seus níveis de expressão
- **Resolução**: Base única

#### 2. `CAGE`
- **Descrição**: Cap Analysis of Gene Expression
- **Uso**: Identificar início de transcrição (promotores)
- **Resolução**: Base única

#### 3. `PROCAP`
- **Descrição**: Precision Run-On sequencing (PRO-cap)
- **Uso**: Mapear início de transcrição com alta precisão
- **Resolução**: Base única

### 🔬 Acessibilidade de Cromatina

#### 4. `ATAC`
- **Descrição**: Assay for Transposase-Accessible Chromatin
- **Uso**: Identificar regiões de cromatina aberta/acessível
- **Resolução**: Base única
- **Aplicação**: Encontrar elementos regulatórios ativos

#### 5. `DNASE`
- **Descrição**: DNase I hypersensitivity
- **Uso**: Identificar regiões de cromatina acessível
- **Resolução**: Base única
- **Similar ao**: ATAC-seq

### ⚛️ Modificações de Cromatina

#### 6. `CHIP_HISTONE`
- **Descrição**: ChIP-seq para marcadores de histonas
- **Inclui**:
  - H3K27AC - Enhancers ativos
  - H3K4ME3 - Promotores ativos
  - H3K27ME3 - Repressão (Polycomb)
  - H3K36ME3 - Corpos gênicos ativos
  - H3K9ME3 - Heterocromatina
  - E outros marcadores
- **Uso**: Identificar estado epigenético da cromatina
- **Resolução**: Base única

#### 7. `CHIP_TF`
- **Descrição**: ChIP-seq para fatores de transcrição
- **Inclui**:
  - CTCF - Insuladores e loops
  - E outros fatores de transcrição
- **Uso**: Identificar sítios de ligação de TFs
- **Resolução**: Base única

### 🧩 Estrutura 3D e Splicing

#### 8. `CONTACT_MAPS`
- **Descrição**: Mapas de contato de cromatina (Hi-C like)
- **Uso**: Predizer interações 3D da cromatina
- **Aplicação**: Entender organização nuclear

#### 9. `SPLICE_JUNCTIONS`
- **Descrição**: Predição de junções de splicing
- **Uso**: Identificar exons e introns
- **Aplicação**: Estudar splicing alternativo

#### 10. `SPLICE_SITES`
- **Descrição**: Predição de sítios de splicing (5' e 3')
- **Uso**: Identificar splice donor/acceptor sites
- **Resolução**: Base única

#### 11. `SPLICE_SITE_USAGE`
- **Descrição**: Predição de uso de sítios de splicing
- **Uso**: Quantificar uso de sítios alternativos
- **Aplicação**: Estudar isoformas

---

## 💡 Exemplos de Uso

### Análise Básica (Expressão e Acessibilidade)
```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ CAGE ATAC
```

### Análise Epigenética Completa
```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs ATAC CHIP_HISTONE CHIP_TF
```

### Análise de Splicing
```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ SPLICE_JUNCTIONS SPLICE_SITES SPLICE_SITE_USAGE
```

### Análise 3D
```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs CONTACT_MAPS CHIP_TF
```

### Análise Completa
```bash
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF DNASE
```

---

## 🔍 Como Verificar Outputs Disponíveis

Execute o script de verificação:
```bash
python check_alphagenome_outputs.py
```

Este script listará todos os outputs que sua instalação do AlphaGenome suporta.

---

## 📋 Comparação com Nomenclatura Antiga

Se você viu documentação com nomes diferentes, aqui está a correspondência:

| Nome Antigo | Nome Correto | Tipo |
|-------------|--------------|------|
| H3K27AC | CHIP_HISTONE | Marcador de histona |
| H3K4ME3 | CHIP_HISTONE | Marcador de histona |
| H3K27ME3 | CHIP_HISTONE | Marcador de histona |
| CTCF | CHIP_TF | Fator de transcrição |

**Nota**: `CHIP_HISTONE` e `CHIP_TF` retornam predições para **múltiplos** marcadores/fatores, não apenas um.

---

## 🎯 Recomendações por Caso de Uso

### 1. Análise de Variantes Regulatórias
```bash
--outputs RNA_SEQ CAGE ATAC CHIP_HISTONE
```
Identifica impacto em expressão, promotores e enhancers.

### 2. Estudos de Expressão Gênica
```bash
--outputs RNA_SEQ CAGE PROCAP
```
Caracteriza expressão e início de transcrição.

### 3. Análise de Elementos Regulatórios
```bash
--outputs ATAC DNASE CHIP_HISTONE CHIP_TF
```
Identifica enhancers, promotores e sítios de ligação.

### 4. Estudos de Splicing
```bash
--outputs RNA_SEQ SPLICE_JUNCTIONS SPLICE_SITES SPLICE_SITE_USAGE
```
Analisa padrões de splicing e isoformas.

### 5. Genômica 3D
```bash
--outputs CONTACT_MAPS CHIP_TF
```
Estuda organização nuclear e loops.

---

## ⚠️ Notas Importantes

1. **Tempo de Processamento**: Quanto mais outputs, mais tempo levará a análise
2. **Tamanho dos Arquivos**: Alguns outputs geram arquivos grandes (especialmente CONTACT_MAPS)
3. **Resolução**: A maioria oferece predições em resolução de base única
4. **Combinações**: Você pode combinar quantos outputs quiser

---

## 🆘 Troubleshooting

### Erro: "Output 'X' não disponível"
**Causa**: Output não existe ou nome incorreto  
**Solução**: Execute `python check_alphagenome_outputs.py` para ver lista correta

### Erro: "Nenhum output válido especificado"
**Causa**: Todos os outputs solicitados são inválidos  
**Solução**: Use outputs da lista acima

### Processo muito lento
**Causa**: Muitos outputs ou sequência muito longa  
**Solução**: Reduza número de outputs ou divida a sequência

---

**Para mais informações**: https://www.alphagenomedocs.com/

*Última atualização: Outubro 2025*

