# Guia: Como Descobrir Tecidos/Células Disponíveis no AlphaGenome

## 🔍 Nova Funcionalidade: `--list-tissues`

O script `build_window_and_predict.py` agora inclui uma funcionalidade para listar todas as ontologias de tecidos/células disponíveis no AlphaGenome!

## 📋 Comandos Disponíveis

### 1. Listar TODOS os tecidos/células

```bash
python3 ~/genomics/build_window_and_predict.py --list-tissues
```

**Saída esperada:**
```
[INFO] Loading tissue metadata from AlphaGenome (this may take a few seconds)...

================================================================================
Available tissues/cells in AlphaGenome (XXX total)
================================================================================

CURIE                     Biosample Name                                Type           
------------------------- --------------------------------------------- ---------------
CL:0000182                hepatocyte                                    primary cell   
CL:0000540                neuron                                        primary cell   
...
UBERON:0000955            brain                                         tissue         
UBERON:0002048            lung                                          tissue         
UBERON:0002107            liver                                         tissue         
...

================================================================================
Usage: --tissue CURIE (e.g., --tissue UBERON:0002107)
================================================================================
```

### 2. Filtrar por nome (RECOMENDADO)

Procurar apenas tecidos/células que contêm "brain":

```bash
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue brain
```

Outros exemplos úteis:
```bash
# Procurar por fígado
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue liver

# Procurar por coração
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue heart

# Procurar por pulmão
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue lung

# Procurar por células T
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue "T cell"

# Procurar por neurônio
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue neuron
```

### 3. Listar tipos de output disponíveis

```bash
python3 ~/genomics/build_window_and_predict.py --list-outputs
```

**Saída:**
```
Available OutputType attributes in AlphaGenome:
  ATAC
  CAGE
  CHIP_HISTONE
  CHIP_TF
  CONTACT_MAPS
  DNASE
  PROCAP
  RNA_SEQ
  SPLICE_JUNCTIONS
  SPLICE_SITES
  SPLICE_SITE_USAGE
```

## 🎯 Usando os CURIEs nas Predições

Depois de encontrar o CURIE desejado, use-o no comando de predição:

```bash
python3 ~/genomics/build_window_and_predict.py \
  --sample HG00096 \
  --gene CYP2B6 \
  --ref-fasta refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --vcf longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --outdir ./alphagenome \
  --predict \
  --outputs "CAGE" \
  --tissue "UBERON:0002107"
```

## 📚 Tipos de Ontologias

O AlphaGenome usa ontologias padronizadas:

| Prefixo | Nome | Descrição | Exemplos |
|---------|------|-----------|----------|
| **UBERON** | Uber-anatomy ontology | Anatomia e tecidos | UBERON:0002107 (liver), UBERON:0000955 (brain) |
| **CL** | Cell Ontology | Tipos celulares | CL:0000182 (hepatocyte), CL:0000540 (neuron) |
| **CLO** | Cell Line Ontology | Linhagens celulares | CLO:0000001 (HeLa) |
| **EFO** | Experimental Factor Ontology | Fatores experimentais | EFO:0000001 |
| **NTR** | New Term Requested | Termos novos solicitados | NTR:XXX |

## 🌐 Navegadores Online de Ontologias

Se preferir navegar visualmente:

### UBERON (Tecidos/Anatomia)
- **URL**: https://www.ebi.ac.uk/ols/ontologies/uberon
- **Uso**: Procure por órgãos e tecidos
- **Exemplos comuns**:
  - Fígado: UBERON:0002107
  - Cérebro: UBERON:0000955
  - Coração: UBERON:0000948
  - Pulmão: UBERON:0002048
  - Rim: UBERON:0002113
  - Pâncreas: UBERON:0001264
  - Baço: UBERON:0002106
  - Sangue: UBERON:0000178

### Cell Ontology (Tipos Celulares)
- **URL**: https://www.ebi.ac.uk/ols/ontologies/cl
- **Uso**: Procure por tipos específicos de células
- **Exemplos comuns**:
  - Hepatócito: CL:0000182
  - Neurônio: CL:0000540
  - Cardiomiócito: CL:0000746
  - Célula T: CL:0000084
  - Célula B: CL:0000236
  - Macrófago: CL:0000235
  - Fibroblasto: CL:0000057

## 💡 Dicas de Uso

### Dica 1: Sempre filtrar por nome primeiro
Evite listar todos os tecidos de uma vez (são centenas!). Use `--filter-tissue`:

```bash
# ❌ Ruim: Lista tudo (centenas de linhas)
python3 build_window_and_predict.py --list-tissues

# ✅ Bom: Lista apenas o que interessa
python3 build_window_and_predict.py --list-tissues --filter-tissue brain
```

### Dica 2: Salvar a lista completa para referência
```bash
python3 ~/genomics/build_window_and_predict.py --list-tissues > tissues_complete.txt
```

### Dica 3: Combinar com grep para busca avançada
```bash
# Procurar tecidos relacionados a sistema nervoso
python3 ~/genomics/build_window_and_predict.py --list-tissues | grep -i nerve

# Procurar células do sistema imune
python3 ~/genomics/build_window_and_predict.py --list-tissues | grep -i "immune\|lymph\|T cell"
```

### Dica 4: Para análises em múltiplos tecidos
Se você precisa comparar múltiplos tecidos, execute o script múltiplas vezes com diferentes `--tissue`:

```bash
# Fígado
python3 build_window_and_predict.py ... --tissue "UBERON:0002107" --outdir ./results_liver

# Cérebro
python3 build_window_and_predict.py ... --tissue "UBERON:0000955" --outdir ./results_brain

# Coração
python3 build_window_and_predict.py ... --tissue "UBERON:0000948" --outdir ./results_heart
```

## ⚠️ Importante: API Key Necessária

A opção `--list-tissues` requer uma API key válida do AlphaGenome porque precisa conectar ao servidor para obter os metadados.

**Formas de fornecer a API key:**

1. Via variável de ambiente (recomendado):
   ```bash
   export ALPHAGENOME_API_KEY="your-key-here"
   python3 build_window_and_predict.py --list-tissues
   ```

2. Via argumento:
   ```bash
   python3 build_window_and_predict.py --list-tissues --api-key "your-key-here"
   ```

## 📊 Estatísticas por Output Type

Número de biosamples e tracks por tipo de output:

| Output Type | Biosamples Únicos | Total de Tracks |
|-------------|-------------------|-----------------|
| RNA_SEQ | 285 | 667 |
| CAGE | 264 | 546 |
| DNASE | 305 | 305 |
| ATAC | 167 | 167 |
| CHIP_HISTONE | 219 | 1116 |
| PROCAP | 6 | 12 |

## 🚀 Workflow Completo

```bash
# 1. Ver outputs disponíveis
python3 build_window_and_predict.py --list-outputs

# 2. Procurar tecido de interesse
python3 build_window_and_predict.py --list-tissues --filter-tissue liver

# 3. Copiar o CURIE desejado (ex: UBERON:0002107)

# 4. Executar predição
python3 build_window_and_predict.py \
  --sample HG00096 \
  --gene CYP2B6 \
  --ref-fasta refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --vcf longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --outdir ./alphagenome \
  --predict \
  --outputs "CAGE" \
  --tissue "UBERON:0002107"

# 5. Analisar resultados
python3 read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/cage.npz
```

## 📝 Referências

- [AlphaGenome Documentation](https://alphafold.com/alphagenome)
- [UBERON Ontology Browser](https://www.ebi.ac.uk/ols/ontologies/uberon)
- [Cell Ontology Browser](https://www.ebi.ac.uk/ols/ontologies/cl)
- [Ontology Lookup Service (OLS)](https://www.ebi.ac.uk/ols/)

