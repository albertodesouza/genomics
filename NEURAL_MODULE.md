# 🧬 Neural Module - Análise de DNA com AlphaGenome AI

## 📖 Visão Geral

O **Neural Module** é uma implementação completa para análise de DNA usando a API do [AlphaGenome](https://github.com/google-deepmind/alphagenome) da Google DeepMind. Este módulo permite predizer características funcionais do DNA através de inteligência artificial, incluindo:

- 🧬 **Expressão Gênica** (RNA-seq, CAGE, PRO-cap)
- 🔬 **Acessibilidade de Cromatina** (ATAC-seq, DNase-seq)
- ⚛️ **Marcadores Epigenéticos** (H3K27AC, H3K4ME3, H3K27ME3, etc.)
- 🔗 **Fatores de Transcrição** (CTCF e outros)
- 🧩 **Estrutura 3D** (Contact Maps)
- ✂️ **Splicing** (Junction sites, site usage)

## 🎯 Principais Características

✅ **11 tipos de análises** suportadas pelo AlphaGenome  
✅ **Visualizações avançadas** (heatmaps, dashboards, comparações)  
✅ **Análise de variantes** com predição de efeitos funcionais  
✅ **Interface de linha de comando** completa e intuitiva  
✅ **Uso programático** como biblioteca Python  
✅ **Integração** com pipelines genômicos existentes  
✅ **Documentação completa** em português e inglês  

---

## 📚 Documentação

### 🚀 Início Rápido
- **[Guia de Instalação](INSTALL_NEURAL.md)** - Como instalar e configurar
- **[Download de Sequências](DOWNLOAD_SEQUENCES.md)** - Como baixar sequências genômicas reais
- **[Guia de Uso](USAGE_NEURAL.md)** - Como executar análises
- **[Interpretando Resultados](RESULTS_NEURAL.md)** - Como interpretar as visualizações

### 📖 Documentação Detalhada
- **[README Completo](NEURAL_MODULE_README.md)** - Documentação técnica completa
- **[Quick Start](NEURAL_QUICKSTART.md)** - Início em 5 minutos
- **[Outputs Disponíveis](OUTPUTS_DISPONIVEIS.md)** - Lista de todos os tipos de análise
- **[Tamanhos Suportados](TAMANHOS_SUPORTADOS.md)** - Restrições de tamanho de sequência

### 🔧 Recursos Avançados
- **[Visualizações Avançadas](VISUALIZACOES_AVANCADAS.md)** - Guia de visualizações
- **[Uso Programático](neural_example.py)** - Exemplos de código Python
- **[Integração com Pipeline](neural_integration.py)** - Ponte com genomes_analyzer

### 🐛 Troubleshooting
- **[Correções Aplicadas](CORRECOES_APLICADAS.md)** - Problemas resolvidos
- **[FAQ](NEURAL_MODULE_README.md#troubleshooting)** - Perguntas frequentes

---

## 💡 Exemplo de Uso

### Análise Básica
```bash
python neural_module.py \
    -i example_sequence.fasta \
    -k YOUR_API_KEY \
    -o results/
```

### Análise de Variante (Anemia Falciforme)
```bash
python neural_module.py \
    -i example_sickle_cell.fasta \
    -k YOUR_API_KEY \
    -o sickle_cell_analysis/ \
    --variant 1024 A T
```

### Análise Completa com Visualizações Avançadas
```bash
python neural_module.py \
    -i gene_region.fasta \
    -k YOUR_API_KEY \
    -o comprehensive_analysis/ \
    --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF \
    --dpi 600 \
    --formats png pdf
```

---

## 🧪 Exemplo Incluído: Anemia Falciforme

O arquivo `example_sequence.fasta` contém a região do gene **HBB** (Beta-globina) que, quando mutado, causa anemia falciforme. Esta é uma das doenças genéticas mais estudadas e serve como excelente exemplo de como uma única mutação pode afetar a função do gene.

**Mutação**: Posição 1024: A→T (GAG→GTG, Glu→Val)

Para analisar:
```bash
python neural_module.py \
    -i example_sequence.fasta \
    -k YOUR_API_KEY \
    -o sickle_cell/ \
    --variant 1024 A T
```

---

## 📊 Saídas Geradas

### Visualizações (Modo Padrão: Avançado)
- **Enhanced tracks** - Múltiplos subplots com metadados
- **Heatmaps** - Comparação de múltiplas tracks
- **Dashboard** - Resumo estatístico completo
- **Comparação multi-output** - Todos outputs em um gráfico
- **Análise de variantes** - 3 painéis (sobreposição, diferença, zoom)

### Relatórios
- **analysis_report.json** - Resumo em JSON
- **Metadados de ontologia** - Informações de tecidos/células

---

## 🗂️ Estrutura de Arquivos

```
genomics/
├── neural_module.py                    # 🌟 Módulo principal
├── neural_example.py                   # 📝 Exemplos de uso
├── neural_integration.py               # 🔗 Integração com pipeline
├── neural_visualizations_advanced.py   # 🎨 Visualizações avançadas
│
├── example_sequence.fasta              # 🧬 Exemplo: Gene HBB (anemia falciforme)
├── neural_config.yaml                  # ⚙️ Configuração
│
├── NEURAL_MODULE.md                    # 📖 Este arquivo
├── INSTALL_NEURAL.md                   # 🚀 Guia de instalação
├── USAGE_NEURAL.md                     # 💡 Guia de uso
├── RESULTS_NEURAL.md                   # 📊 Interpretando resultados
│
└── (outros arquivos de documentação)
```

---

## 🔗 Links Úteis

- **Documentação do AlphaGenome**: https://www.alphagenomedocs.com/
- **GitHub do AlphaGenome**: https://github.com/google-deepmind/alphagenome
- **Obter API Key**: https://www.alphagenomedocs.com/
- **Paper**: Avsec et al. 2025 - "AlphaGenome: advancing regulatory variant effect prediction"

---

## 🤝 Contribuindo

O Neural Module é parte do projeto de análise genômica. Para contribuir:
1. Reporte bugs via issues
2. Sugira melhorias
3. Compartilhe casos de uso

---

## 📄 Licença

Compatível com Apache 2.0 (licença do AlphaGenome)  
Uso gratuito para pesquisa não comercial

---

## 📞 Suporte

- **AlphaGenome**: alphagenome@google.com
- **Documentação**: Veja os guias listados acima
- **Issues**: Abra uma issue no repositório

---

## ✨ Status do Projeto

🟢 **Estável e Pronto para Uso**

- ✅ Totalmente testado
- ✅ Documentação completa
- ✅ Exemplos funcionando
- ✅ Visualizações profissionais
- ✅ Integração com pipeline

---

**Desenvolvido com ❤️ para análise genômica avançada**

*Última atualização: Outubro 2025*

