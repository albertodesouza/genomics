# 🚀 Guia de Instalação - Neural Module

## 📋 Pré-requisitos

### Sistema Operacional
- ✅ Linux (Ubuntu, Debian, CentOS, etc.)
- ✅ macOS
- ⚠️ Windows (via WSL2 recomendado)

### Software Necessário
- **Python**: 3.7 ou superior
- **Conda**: Miniconda ou Anaconda (recomendado)
- **Git**: Para clonar repositórios

---

## 📦 Instalação Passo a Passo

### Passo 1: Verificar Requisitos

Execute o script de verificação:
```bash
bash check_neural_requirements.sh
```

Este script verifica:
- ✓ Python e versão
- ✓ Módulos necessários (rich, matplotlib, numpy)
- ✓ AlphaGenome instalado
- ✓ Conectividade com internet

---

### Passo 2: Instalar AlphaGenome

#### Opção A: Script Automático (Recomendado)
```bash
bash install_alphagenome.sh
```

#### Opção B: Manual
```bash
git clone https://github.com/google-deepmind/alphagenome.git
pip install ./alphagenome
```

---

### Passo 3: Instalar Dependências Python

```bash
# Dependências obrigatórias
pip install rich matplotlib numpy pandas

# Dependências opcionais (para neural_integration.py)
pip install pyyaml
```

---

### Passo 4: Obter API Key

1. Acesse: **https://www.alphagenomedocs.com/**
2. Crie uma conta (gratuito para uso não comercial)
3. Solicite sua API key
4. Guarde a chave em local seguro

**⚠️ IMPORTANTE**: Nunca compartilhe sua API key publicamente!

---

### Passo 5: Testar Instalação

```bash
# Teste básico
python neural_module.py --help

# Teste completo (requer API key)
bash test_neural_module.sh YOUR_API_KEY
```

Se tudo funcionar, você verá:
```
✓ 1/1 sequências processadas com sucesso
```

---

## 🐍 Instalação com Conda (Recomendado)

### Criar Ambiente

```bash
# Criar ambiente
conda create -n neural python=3.10

# Ativar
conda activate neural

# Instalar dependências
conda install -c conda-forge matplotlib numpy pandas rich pyyaml

# Instalar AlphaGenome
bash install_alphagenome.sh
```

### Uso Futuro
```bash
conda activate neural
python neural_module.py -i seq.fasta -k API_KEY -o results/
```

---

## 🔧 Ferramentas Opcionais

Para usar `neural_integration.py` (extração de sequências de VCF/BED):

```bash
# Via conda
conda install -c bioconda bcftools samtools bedtools

# Via apt (Ubuntu/Debian)
sudo apt install bcftools samtools bedtools

# Via homebrew (macOS)
brew install bcftools samtools bedtools
```

---

## ✅ Verificação Final

Execute este checklist:

```bash
# 1. Python disponível?
python3 --version
# Deve mostrar: Python 3.x.x

# 2. Módulos instalados?
python3 -c "import rich, matplotlib, numpy; print('✓ OK')"
# Deve mostrar: ✓ OK

# 3. AlphaGenome instalado?
python3 -c "from alphagenome.models import dna_client; print('✓ OK')"
# Deve mostrar: ✓ OK

# 4. Neural module funcionando?
python3 neural_module.py --help
# Deve mostrar a ajuda

# 5. Exemplo de sequência existe?
ls -lh example_sequence.fasta
# Deve mostrar o arquivo
```

---

## 🐛 Problemas Comuns

### Erro: "No module named 'rich'"
```bash
pip install rich
```

### Erro: "No module named 'alphagenome'"
```bash
bash install_alphagenome.sh
```

### Erro: "ModuleNotFoundError: No module named 'matplotlib'"
```bash
pip install matplotlib
```

### Erro: Conexão com API falhou
- Verifique sua conexão com internet
- Confirme que a API key está correta
- Teste: `ping google.com`

---

## 📊 Espaço em Disco

Requisitos de espaço:
- **AlphaGenome**: ~500 MB
- **Dependências Python**: ~200 MB
- **Neural Module**: ~5 MB
- **Resultados** (por análise): ~5-50 MB

**Total recomendado**: 2 GB livres

---

## 🔄 Atualizações

### Atualizar AlphaGenome
```bash
pip install --upgrade alphagenome
```

### Atualizar Dependências
```bash
pip install --upgrade rich matplotlib numpy pandas
```

---

## 📚 Próximos Passos

Após instalação bem-sucedida:

1. ✅ **[Guia de Uso](USAGE_NEURAL.md)** - Aprenda a executar análises
2. ✅ **[Quick Start](NEURAL_QUICKSTART.md)** - Exemplo rápido
3. ✅ **[Exemplo de Anemia Falciforme](NEURAL_MODULE.md#exemplo-incluído-anemia-falciforme)** - Caso de uso real

---

## 💡 Dicas

1. **Use ambiente virtual** para evitar conflitos
2. **Guarde API key** em local seguro (não no código!)
3. **Teste primeiro** com sequência de exemplo
4. **Leia a documentação** do AlphaGenome

---

## 🆘 Suporte

Se encontrar problemas:

1. Execute: `bash check_neural_requirements.sh`
2. Consulte: [CORRECOES_APLICADAS.md](CORRECOES_APLICADAS.md)
3. Veja: [FAQ no README](NEURAL_MODULE_README.md#troubleshooting)

---

**Instalação concluída? Vá para**: [Guia de Uso →](USAGE_NEURAL.md)

*Última atualização: Outubro 2025*

