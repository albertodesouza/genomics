# 📁 Estrutura do Módulo build_non_longevous_dataset

Este diretório contém o módulo **Non-Longevous Dataset Builder** completamente organizado.

## 📂 Estrutura de Arquivos

```
build_non_longevous_dataset/
├── build_non_longevous_dataset.py    # Programa principal
├── build_window_and_predict.py       # Script para extrair janelas genômicas e predições
├── README.md                         # Documentação completa
├── QUICKSTART.md                     # Guia de início rápido
├── IMPLEMENTACAO.md                  # Detalhes técnicos da implementação
├── ESTRUTURA.md                      # Este arquivo
├── configs/
│   └── default.yaml                  # Configuração padrão
├── docs/
│   ├── ALPHAGENOME_PREDICTIONS.md    # Guia de predições AlphaGenome
│   └── ALPHAGENOME_TISSUES.md        # Guia de tecidos/células AlphaGenome
└── scripts/
    └── test.sh                       # Script de teste
```

## 🚀 Como Usar

### Do diretório do módulo:
```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

### Da raiz do projeto:
```bash
python3 build_non_longevous_dataset/build_non_longevous_dataset.py \
  --config build_non_longevous_dataset/configs/default.yaml
```

### Usando o script de teste:
```bash
cd build_non_longevous_dataset/scripts
bash test.sh
```

## 📝 Caminhos Relativos

Os caminhos no arquivo `configs/default.yaml` são relativos ao diretório `configs/`:
- `../../doc/arquivo.csv` → `/caminho/para/genomics/doc/arquivo.csv`
- `../../refs/genoma.fa` → `/caminho/para/genomics/refs/genoma.fa`

## 📚 Documentação

- **README.md**: Documentação completa do módulo
- **QUICKSTART.md**: Guia rápido para começar
- **IMPLEMENTACAO.md**: Detalhes técnicos da implementação

## ✅ Testado e Funcionando

✓ Execução do diretório do módulo
✓ Execução da raiz do projeto
✓ Script de teste funcional
✓ Resolução correta de caminhos relativos
✓ Integração com build_window_and_predict.py (incluído no módulo)

## 🔗 build_window_and_predict.py

Este módulo inclui o `build_window_and_predict.py`, que é responsável por:
- Extrair janelas genômicas de 1 Mb ao redor de genes específicos
- Aplicar variantes de amostras do 1000 Genomes
- Gerar sequências consenso por haplótipo
- Executar predições AlphaGenome (opcional)

📚 Documentação adicional:
- [Guia de Predições AlphaGenome](docs/ALPHAGENOME_PREDICTIONS.md)
- [Guia de Tecidos/Células](docs/ALPHAGENOME_TISSUES.md)
