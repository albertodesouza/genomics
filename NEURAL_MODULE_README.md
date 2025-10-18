# Neural Module - Análise de DNA com AlphaGenome

## 📋 Descrição

O `neural_module.py` é uma ferramenta Python que integra a API do [AlphaGenome](https://github.com/google-deepmind/alphagenome) da Google DeepMind para realizar análises avançadas de sequências de DNA.

## 🚀 Funcionalidades

- **Análise de Sequências**: Predição de múltiplos aspectos funcionais de sequências de DNA
- **Predição de Expressão Gênica**: RNA-seq, CAGE
- **Características de Cromatina**: ATAC-seq, H3K27AC, H3K4ME3, H3K27ME3, H3K36ME3, H3K9ME3, CTCF
- **Análise de Variantes**: Predição do efeito de variantes SNP
- **Visualizações**: Geração automática de gráficos em múltiplos formatos (PNG, PDF, SVG)
- **Suporte a FASTA**: Leitura de arquivos FASTA padrão

## 📦 Instalação

### 1. Instalar o AlphaGenome

Execute o script de instalação fornecido:

```bash
bash install_alphagenome.sh
```

Ou manualmente:

```bash
git clone https://github.com/google-deepmind/alphagenome.git
pip install ./alphagenome
```

### 2. Instalar Dependências Adicionais

```bash
pip install rich matplotlib
```

### 3. Obter API Key

Acesse [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/) e solicite sua API key gratuita para uso não comercial.

## 🎯 Uso

### Exemplo Básico

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/
```

### Análise com Outputs Específicos

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --outputs RNA_SEQ ATAC H3K27AC
```

### Análise de Variante

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --variant 1000 A C
```

Este comando analisa o efeito de uma variante A→C na posição 1000 (relativa ao início da sequência).

### Múltiplos Formatos de Saída

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --formats png pdf svg --dpi 600
```

### Apenas Análise (Sem Gráficos)

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --no-plots
```

## 📊 Tipos de Output Disponíveis

| Output | Descrição |
|--------|-----------|
| `RNA_SEQ` | Predição de expressão gênica via RNA-seq |
| `CAGE` | Cap Analysis of Gene Expression |
| `ATAC` | Acessibilidade de cromatina (ATAC-seq) |
| `H3K27AC` | Marcador de elementos regulatórios ativos |
| `H3K4ME3` | Marcador de promotores ativos |
| `H3K27ME3` | Marcador de repressão gênica |
| `H3K36ME3` | Marcador de corpos gênicos ativos |
| `H3K9ME3` | Marcador de heterocromatina |
| `CTCF` | Fator de ligação de insulador |

## 📝 Formato de Entrada (FASTA)

O arquivo FASTA deve seguir o formato padrão:

```
>sequence_id_1 description
ATCGATCGATCGATCG...
>sequence_id_2 description
GCTAGCTAGCTAGCTA...
```

### Requisitos:

- Sequências com 100 bp a 1.000.000 bp (1 Mbp)
- Caracteres válidos: A, C, G, T, N (e códigos de ambiguidade IUPAC)

## 📁 Estrutura de Saída

```
results/
├── sequence_id_1_RNA_SEQ.png
├── sequence_id_1_ATAC.png
├── sequence_id_1_H3K27AC.png
├── sequence_id_2_RNA_SEQ.png
├── ...
└── analysis_report.json
```

### Relatório JSON

O arquivo `analysis_report.json` contém:

```json
{
  "timestamp": "2025-10-16T10:30:00",
  "total_sequences": 2,
  "successful_analyses": 2,
  "sequences": [
    {
      "id": "sequence_id_1",
      "length": 50000,
      "status": "success",
      "outputs": ["RNA_SEQ", "ATAC", "H3K27AC"]
    }
  ]
}
```

## ⚙️ Opções de Linha de Comando

### Obrigatórias

- `-i, --input`: Arquivo FASTA de entrada
- `-k, --api-key`: Chave da API do AlphaGenome
- `-o, --output`: Diretório de saída

### Opcionais

- `--outputs`: Tipos de output desejados (padrão: RNA_SEQ CAGE ATAC H3K27AC H3K4ME3)
- `--chromosome`: Cromossomo de referência (padrão: chr1)
- `--start`: Posição inicial de referência (padrão: 1000000)
- `--variant POS REF ALT`: Analisar variante na posição POS com bases REF>ALT
- `--formats`: Formatos de gráficos (png, pdf, svg) (padrão: png)
- `--dpi`: Resolução dos gráficos (padrão: 300)
- `--no-plots`: Não gerar gráficos (apenas análise)

## 💡 Exemplos Avançados

### 1. Análise Completa de uma Região Genômica

```bash
python neural_module.py \
    -i chr1_region.fasta \
    -k YOUR_API_KEY \
    -o chr1_analysis/ \
    --chromosome chr1 \
    --start 1000000 \
    --outputs RNA_SEQ CAGE ATAC H3K27AC H3K4ME3 H3K27ME3 CTCF \
    --formats png pdf \
    --dpi 600
```

### 2. Análise Rápida de Múltiplas Sequências

```bash
python neural_module.py \
    -i multiple_sequences.fasta \
    -k YOUR_API_KEY \
    -o quick_analysis/ \
    --outputs RNA_SEQ ATAC \
    --formats png \
    --dpi 150
```

### 3. Análise de Variante Patogênica

```bash
python neural_module.py \
    -i disease_gene.fasta \
    -k YOUR_API_KEY \
    -o variant_effect/ \
    --variant 5000 G T \
    --formats png pdf svg
```

## 🔧 Troubleshooting

### Erro: "AlphaGenome não está instalado"

```bash
pip install git+https://github.com/google-deepmind/alphagenome.git
```

### Erro: "API key inválida"

Verifique se sua API key está correta e ativa em [alphagenomedocs.com](https://www.alphagenomedocs.com/)

### Erro: "Sequência muito longa"

AlphaGenome suporta até 1 Mbp. Divida sequências maiores em pedaços menores.

### Problema com Matplotlib

```bash
pip install --upgrade matplotlib seaborn
```

## 📚 Recursos

- **Documentação do AlphaGenome**: [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/)
- **Repositório GitHub**: [https://github.com/google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome)
- **Paper**: Avsec et al. 2025 - "AlphaGenome: advancing regulatory variant effect prediction"

## 🔬 Casos de Uso

1. **Análise de Variantes Regulatórias**: Identificar impacto de SNPs em regiões regulatórias
2. **Predição de Elementos Funcionais**: Identificar promotores, enhancers, insuladores
3. **Estudos de Expressão Gênica**: Predizer níveis de expressão em diferentes tecidos
4. **Análise de Cromatina**: Estudar acessibilidade e modificações de histonas
5. **Genômica Funcional**: Caracterizar sequências de DNA desconhecidas

## ⚠️ Limitações

- Requer conexão com internet (API online)
- Taxa de queries limitada (verificar termos de uso)
- Uso gratuito apenas para pesquisa não comercial
- Sequências de 100 bp a 1 Mbp
- Tempo de processamento varia com tamanho da sequência

## 📄 Licença

Este módulo é distribuído sob licença Apache 2.0, compatível com o AlphaGenome.

## 🤝 Contribuições

Contribuições são bem-vindas! Por favor, abra issues ou pull requests no repositório do projeto.

## 📧 Suporte

Para questões sobre o AlphaGenome: alphagenome@google.com
Para questões sobre este módulo: abra uma issue no repositório

---

**Desenvolvido para integração com o genomes_analyzer.py**

