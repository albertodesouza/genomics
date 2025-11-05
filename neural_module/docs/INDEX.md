# 📚 Neural Module - Índice Completo de Arquivos

## 🎯 Visão Geral

O **Neural Module** é uma implementação completa para análise de DNA usando a API do AlphaGenome (Google DeepMind). Este documento serve como índice de todos os arquivos criados e suas funções.

---

## 📁 Arquivos Principais

### 1. `neural_module.py` ⭐
**Descrição**: Módulo principal de análise de DNA com AlphaGenome

**Funcionalidades**:
- Parse de arquivos FASTA
- Validação de sequências de DNA
- Integração com API do AlphaGenome
- Predição de características funcionais (RNA-seq, ATAC, marcadores de histonas)
- Análise de efeitos de variantes
- Geração automática de visualizações
- Relatórios em JSON

**Uso**:
```bash
python neural_module.py -i input.fasta -k API_KEY -o output_dir/
```

**Características**:
- ✅ Interface de linha de comando completa
- ✅ Suporte a múltiplos formatos de saída (PNG, PDF, SVG)
- ✅ Validação robusta de entrada
- ✅ Logging rico com Rich library
- ✅ ~650 linhas de código bem documentado

---

### 2. `neural_example.py`
**Descrição**: Exemplos de uso programático do neural_module

**Exemplos Incluídos**:
1. Análise básica de sequência
2. Análise de arquivo FASTA
3. Análise de variante
4. Análise múltipla (batch)
5. Configuração customizada
6. Integração com genomes_analyzer
7. Processamento em lote de diretórios

**Uso**:
```bash
# Todos os exemplos
python neural_example.py -k API_KEY

# Exemplo específico
python neural_example.py -k API_KEY -e 3
```

---

### 3. `neural_integration.py`
**Descrição**: Ponte entre genomes_analyzer.py e neural_module.py

**Funcionalidades**:
- Extração de sequências de arquivos VCF
- Extração de sequências de arquivos BED
- Extração de sequências de genes específicos (via GTF)
- Análise integrada completa (VCF → FASTA → AlphaGenome)
- Correlação de resultados

**Uso**:
```bash
# Análise integrada
python neural_integration.py \
    --integrated \
    --vcf variants.vcf \
    --ref genome.fa \
    --api-key API_KEY \
    --output results/

# Extrair sequências de genes
python neural_integration.py \
    --extract-genes \
    --genes BRCA1 TP53 \
    --gtf genes.gtf \
    --ref genome.fa \
    --output genes.fasta
```

**Requer**: bcftools, samtools, bedtools (opcionais)

---

## 🛠️ Scripts de Instalação e Configuração

### 4. `install_alphagenome.sh`
**Descrição**: Script automático para instalação do AlphaGenome

**O que faz**:
- Clona repositório do AlphaGenome
- Instala via pip
- Limpa arquivos temporários
- Verifica ambiente conda (se disponível)

**Uso**:
```bash
bash install_alphagenome.sh
```

---

### 5. `check_neural_requirements.sh`
**Descrição**: Script de diagnóstico completo

**Verificações**:
- ✅ Arquivos do projeto
- ✅ Python e versão
- ✅ Módulos Python (rich, matplotlib, numpy, pandas)
- ✅ AlphaGenome e submódulos
- ✅ Ferramentas opcionais (bcftools, samtools, bedtools)
- ✅ Conectividade com internet
- ✅ Ambiente conda

**Uso**:
```bash
bash check_neural_requirements.sh
```

**Saída**: Relatório detalhado com ✅/❌/⚠️

---

## 🧪 Scripts de Teste

### 6. `test_neural_module.sh`
**Descrição**: Suite de testes automatizados

**Testes Incluídos**:
1. Análise básica com RNA_SEQ
2. Múltiplos outputs
3. Múltiplos formatos de saída
4. Análise sem gráficos

**Uso**:
```bash
bash test_neural_module.sh YOUR_API_KEY
```

**Saída**: Diretório `neural_test_TIMESTAMP/` com resultados

---

### 7. `demo_neural_module.sh`
**Descrição**: Demonstração interativa e guia de uso

**Conteúdo**:
- Pré-requisitos
- Exemplos de uso
- Outputs disponíveis
- Casos de uso
- Opções avançadas
- Recursos e suporte

**Uso**:
```bash
bash demo_neural_module.sh
# ou
cat demo_neural_module.sh
```

---

## 📖 Documentação

### 8. `NEURAL_MODULE_README.md` 📘
**Descrição**: Documentação completa e detalhada

**Seções**:
- Descrição e funcionalidades
- Instalação passo a passo
- Guia de uso completo
- Tipos de output disponíveis
- Formato de entrada (FASTA)
- Estrutura de saída
- Opções de linha de comando
- Exemplos avançados
- Troubleshooting
- Recursos e links
- Casos de uso
- Limitações
- Licença

**Páginas**: ~15 páginas de documentação

---

### 9. `NEURAL_QUICKSTART.md` 🚀
**Descrição**: Guia de início rápido (5 minutos)

**Conteúdo**:
- Início em 5 minutos
- Comandos úteis
- Integração com genomes_analyzer
- Exemplos práticos
- Troubleshooting rápido
- Recursos adicionais

**Para**: Usuários que querem começar rapidamente

---

### 10. `NEURAL_MODULE_INDEX.md` 📚
**Descrição**: Este arquivo - índice completo

---

## ⚙️ Arquivos de Configuração

### 11. `neural_config.yaml`
**Descrição**: Arquivo de configuração YAML

**Configurações**:
- API key do AlphaGenome
- Outputs padrão
- Parâmetros de visualização
- Limites de sequência
- Contexto genômico
- Processamento

**Uso**: Pode ser usado para armazenar configurações persistentes

---

## 📄 Arquivos de Exemplo

### 12. `example_sequence.fasta`
**Descrição**: Arquivo FASTA de exemplo para testes

**Conteúdo**:
- 2 sequências de teste
- test_sequence_1: 600 bp
- test_sequence_2: 600 bp

**Uso**: Para testar o neural_module sem precisar de dados reais

---

## 📊 Estrutura de Diretórios

```
genomics/
│
├── 🐍 Módulos Python
│   ├── neural_module.py          (Principal)
│   ├── neural_example.py         (Exemplos)
│   └── neural_integration.py     (Integração)
│
├── 🔧 Scripts Shell
│   ├── install_alphagenome.sh
│   ├── check_neural_requirements.sh
│   ├── test_neural_module.sh
│   └── demo_neural_module.sh
│
├── 📖 Documentação
│   ├── NEURAL_MODULE_README.md   (Completa)
│   ├── NEURAL_QUICKSTART.md      (Rápida)
│   └── NEURAL_MODULE_INDEX.md    (Este arquivo)
│
├── ⚙️ Configuração
│   └── neural_config.yaml
│
└── 📄 Exemplos
    └── example_sequence.fasta
```

---

## 🎯 Casos de Uso

### Para Iniciantes
1. Ler `NEURAL_QUICKSTART.md`
2. Executar `check_neural_requirements.sh`
3. Executar `bash install_alphagenome.sh`
4. Testar com `example_sequence.fasta`

### Para Uso Direto
```bash
python neural_module.py -i seu_arquivo.fasta -k API_KEY -o results/
```

### Para Integração com Pipeline
```bash
python neural_integration.py --integrated \
    --vcf resultados/variants.vcf \
    --ref reference.fa \
    --api-key API_KEY \
    --output integrated/
```

### Para Desenvolvimento
- Consultar `neural_example.py` para uso programático
- Importar classes e funções de `neural_module.py`
- Estender funcionalidades conforme necessário

---

## 📈 Estatísticas

| Arquivo | Linhas | Tipo | Status |
|---------|--------|------|--------|
| neural_module.py | ~650 | Python | ✅ Completo |
| neural_example.py | ~450 | Python | ✅ Completo |
| neural_integration.py | ~550 | Python | ✅ Completo |
| install_alphagenome.sh | ~45 | Shell | ✅ Completo |
| check_neural_requirements.sh | ~250 | Shell | ✅ Completo |
| test_neural_module.sh | ~90 | Shell | ✅ Completo |
| demo_neural_module.sh | ~200 | Shell | ✅ Completo |
| NEURAL_MODULE_README.md | ~550 | Markdown | ✅ Completo |
| NEURAL_QUICKSTART.md | ~250 | Markdown | ✅ Completo |
| neural_config.yaml | ~80 | YAML | ✅ Completo |
| example_sequence.fasta | ~20 | FASTA | ✅ Completo |

**Total**: ~3.135 linhas de código e documentação

---

## 🚀 Início Rápido

### Método 1: Teste Rápido (Recomendado)
```bash
# 1. Verificar requisitos
bash check_neural_requirements.sh

# 2. Se necessário, instalar AlphaGenome
bash install_alphagenome.sh

# 3. Testar
python neural_module.py \
    -i example_sequence.fasta \
    -k YOUR_API_KEY \
    -o test_results/
```

### Método 2: Teste Completo
```bash
bash test_neural_module.sh YOUR_API_KEY
```

### Método 3: Exemplos Programáticos
```bash
python neural_example.py -k YOUR_API_KEY
```

---

## 🔗 Links Úteis

- **AlphaGenome**: https://github.com/google-deepmind/alphagenome
- **Documentação API**: https://www.alphagenomedocs.com/
- **Obter API Key**: https://www.alphagenomedocs.com/
- **Paper**: Avsec et al. 2025 - "AlphaGenome: advancing regulatory variant effect prediction"

---

## 📞 Suporte

### Para Questões sobre AlphaGenome
- Email: alphagenome@google.com
- Documentação: https://www.alphagenomedocs.com/

### Para Questões sobre Neural Module
- Consulte a documentação neste repositório
- Abra issues no repositório (se aplicável)

---

## ✅ Checklist de Implementação

- [x] Módulo principal (neural_module.py)
- [x] Exemplos de uso (neural_example.py)
- [x] Integração com pipeline (neural_integration.py)
- [x] Script de instalação
- [x] Script de diagnóstico
- [x] Suite de testes
- [x] Demonstração
- [x] Documentação completa
- [x] Guia rápido
- [x] Arquivo de configuração
- [x] Sequências de exemplo
- [x] Este índice

**Status**: ✅ 100% Completo

---

## 🎓 Próximos Passos Sugeridos

1. **Testes**: Execute a suite de testes para verificar funcionamento
2. **Personalização**: Ajuste `neural_config.yaml` conforme suas necessidades
3. **Integração**: Use `neural_integration.py` para conectar ao seu pipeline
4. **Exploração**: Teste diferentes tipos de outputs e visualizações
5. **Desenvolvimento**: Estenda funcionalidades conforme necessário

---

**Desenvolvido para análise genômica avançada com AlphaGenome** 🧬

*Última atualização: Outubro 2025*

