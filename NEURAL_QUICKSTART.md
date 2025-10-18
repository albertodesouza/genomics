# 🚀 Neural Module - Guia de Início Rápido

## ⚡ Início em 5 Minutos

### 1. Verificar Requisitos

```bash
bash check_neural_requirements.sh
```

### 2. Instalar AlphaGenome

```bash
bash install_alphagenome.sh
```

### 3. Obter API Key

Acesse: [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/)

### 4. Executar Análise

```bash
python neural_module.py \
    -i example_sequence.fasta \
    -k YOUR_API_KEY \
    -o results/
```

---

## 📋 Comandos Úteis

### Análise Básica

```bash
# Análise completa com outputs padrão
python neural_module.py -i sequence.fasta -k API_KEY -o results/
```

### Análise Customizada

```bash
# Escolher outputs específicos
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ ATAC H3K27AC

# Alta resolução
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --dpi 600 \
    --formats png pdf svg
```

### Análise de Variante

```bash
# Analisar efeito de variante A→C na posição 1000
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --variant 1000 A C
```

---

## 🔗 Integração com genomes_analyzer.py

### Fluxo de Trabalho Completo

```bash
# 1. Executar pipeline genômico
python genomes_analyzer.py config.yaml

# 2. Extrair sequências de interesse
python neural_integration.py \
    --extract-vcf \
    --vcf results/variants.vcf \
    --ref reference.fa \
    --output extracted.fasta

# 3. Analisar com AlphaGenome
python neural_module.py \
    -i extracted.fasta \
    -k API_KEY \
    -o neural_results/

# 4. OU fazer tudo de uma vez:
python neural_integration.py \
    --integrated \
    --vcf results/variants.vcf \
    --ref reference.fa \
    --api-key API_KEY \
    --output integrated_results/
```

### Extrair Genes Específicos

```bash
python neural_integration.py \
    --extract-genes \
    --genes BRCA1 TP53 EGFR \
    --gtf annotations.gtf \
    --ref reference.fa \
    --output genes.fasta
```

---

## 📊 Tipos de Output

| Output | Descrição | Uso Típico |
|--------|-----------|------------|
| `RNA_SEQ` | Expressão gênica | Identificar genes ativos |
| `CAGE` | Cap analysis | Identificar promotores |
| `ATAC` | Acessibilidade | Regiões regulatórias abertas |
| `H3K27AC` | Marcador de histona | Enhancers ativos |
| `H3K4ME3` | Marcador de histona | Promotores ativos |
| `H3K27ME3` | Marcador de histona | Repressão Polycomb |
| `CTCF` | Fator de transcrição | Isoladores/loops |

---

## 💡 Exemplos de Uso

### Exemplo 1: Pesquisa de Variantes Patogênicas

```bash
# Analisar variante em gene de doença
python neural_module.py \
    -i disease_gene.fasta \
    -k API_KEY \
    -o pathogenic_analysis/ \
    --variant 5234 G A \
    --outputs RNA_SEQ CAGE ATAC H3K27AC \
    --formats png pdf
```

### Exemplo 2: Caracterização de Região Genômica

```bash
# Análise completa de região
python neural_module.py \
    -i genomic_region.fasta \
    -k API_KEY \
    -o region_analysis/ \
    --outputs RNA_SEQ CAGE ATAC H3K27AC H3K4ME3 H3K27ME3 CTCF \
    --dpi 600
```

### Exemplo 3: Comparação de Múltiplas Variantes

```bash
# Criar FASTAs com diferentes variantes
# Depois analisar cada uma
for variant in variant1.fasta variant2.fasta variant3.fasta; do
    python neural_module.py \
        -i $variant \
        -k API_KEY \
        -o comparison/$(basename $variant .fasta)/
done
```

---

## 🧪 Testes

### Teste Rápido

```bash
# Teste com sequência de exemplo
python neural_module.py \
    -i example_sequence.fasta \
    -k API_KEY \
    -o test_results/
```

### Suite de Testes Completa

```bash
bash test_neural_module.sh YOUR_API_KEY
```

### Exemplos Programáticos

```bash
# Executar todos os exemplos
python neural_example.py -k API_KEY

# Executar exemplo específico
python neural_example.py -k API_KEY -e 3
```

---

## 📁 Estrutura de Arquivos

```
genomics/
├── neural_module.py              # Módulo principal
├── neural_example.py             # Exemplos de uso
├── neural_integration.py         # Integração com genomes_analyzer
├── neural_config.yaml            # Configuração
├── install_alphagenome.sh        # Script de instalação
├── check_neural_requirements.sh  # Verificação de requisitos
├── test_neural_module.sh         # Suite de testes
├── demo_neural_module.sh         # Demonstração
├── example_sequence.fasta        # Sequência de exemplo
├── NEURAL_MODULE_README.md       # Documentação completa
└── NEURAL_QUICKSTART.md          # Este guia
```

---

## ❓ Troubleshooting

### Problema: "AlphaGenome não está instalado"

```bash
bash install_alphagenome.sh
```

### Problema: "API key inválida"

- Verifique se copiou a key corretamente
- Certifique-se de que está ativa em alphagenomedocs.com

### Problema: "Sequência muito longa"

- AlphaGenome suporta até 1 Mbp
- Divida sequências maiores em pedaços

### Problema: Gráficos não são gerados

```bash
pip install --upgrade matplotlib seaborn
```

### Problema: ImportError

```bash
# Reinstalar dependências
pip install --force-reinstall rich matplotlib numpy
```

---

## 🔍 Verificar Instalação

```bash
# Verificar neural_module
python -c "from neural_module import AlphaGenomeAnalyzer; print('✓ OK')"

# Verificar AlphaGenome
python -c "from alphagenome.models import dna_client; print('✓ OK')"

# Verificar dependências
python -c "import rich, matplotlib, numpy; print('✓ OK')"
```

---

## 📚 Recursos Adicionais

- **Documentação Completa**: `NEURAL_MODULE_README.md`
- **Demonstração**: `bash demo_neural_module.sh`
- **API AlphaGenome**: [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/)
- **GitHub AlphaGenome**: [https://github.com/google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome)

---

## 💬 Suporte

### Perguntas sobre AlphaGenome
- Email: alphagenome@google.com
- Documentação: https://www.alphagenomedocs.com/

### Perguntas sobre neural_module
- Abra uma issue no repositório
- Consulte a documentação completa

---

## ⚠️ Notas Importantes

1. **Uso Não Comercial**: API gratuita apenas para pesquisa
2. **Taxa Limitada**: Adequado para ~1000s de predições
3. **Internet Necessária**: API requer conexão online
4. **Segurança**: Nunca compartilhe sua API key

---

## 🎯 Próximos Passos

Após executar o teste básico:

1. ✅ Explorar diferentes outputs disponíveis
2. ✅ Testar análise de variantes
3. ✅ Integrar com seu pipeline existente
4. ✅ Experimentar com suas próprias sequências
5. ✅ Usar exemplos programáticos para workflows customizados

**Boa análise! 🧬🚀**

