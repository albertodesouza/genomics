# ⚠️ Tamanhos de Sequências Suportados pelo AlphaGenome

## 🎯 Informação Importante

O AlphaGenome **só aceita sequências com tamanhos específicos**:

| Tamanho | Bytes | Uso Recomendado |
|---------|-------|-----------------|
| **2,048 bp** | 2 KB | Pequenas regiões, testes rápidos |
| **16,384 bp** | 16 KB | Genes individuais |
| **131,072 bp** | 128 KB | Múltiplos genes, clusters |
| **524,288 bp** | 512 KB | Grandes regiões genômicas |
| **1,048,576 bp** | 1 MB | Regiões cromossômicas extensas |

## ❌ O Que Acontece com Outros Tamanhos?

Se você tentar usar uma sequência de 600 bp, 5000 bp, ou qualquer outro tamanho não listado, receberá o erro:

```
ValueError: Sequence length 600 not supported by the model. 
Supported lengths: [2048, 16384, 131072, 524288, 1048576]
```

## ✅ Como Resolver

### Opção 1: Pad (Adicionar Bases)

Se sua sequência é **menor** que o tamanho mínimo (2048 bp), adicione bases N:

```python
sequence = "ATCGATCG..."  # 600 bp
padded = sequence + "N" * (2048 - len(sequence))  # Agora tem 2048 bp
```

### Opção 2: Truncar

Se sua sequência é maior mas não é um tamanho suportado, truncar para o próximo tamanho menor:

```python
sequence = "ATCGATCG..." * 1000  # 5000 bp
truncated = sequence[:2048]  # Usa primeiros 2048 bp
```

### Opção 3: Resize para Próximo Tamanho Suportado

```python
def resize_to_supported(seq):
    supported = [2048, 16384, 131072, 524288, 1048576]
    
    # Encontrar o tamanho suportado mais próximo
    for size in supported:
        if len(seq) <= size:
            if len(seq) < size:
                # Adicionar Ns
                return seq + "N" * (size - len(seq))
            return seq
    
    # Se maior que 1MB, truncar
    return seq[:1048576]
```

## 📝 Exemplos Práticos

### Exemplo 1: Sequência Pequena (600 bp → 2048 bp)

```bash
# Criar sequência com padding em Python
python3 << 'EOF'
seq = "ATCG" * 150  # 600 bp
padded = seq + "N" * (2048 - len(seq))
with open("seq_2kb.fasta", "w") as f:
    f.write(">my_sequence_2kb\n")
    f.write(padded + "\n")
EOF

# Analisar
python neural_module.py -i seq_2kb.fasta -k API_KEY -o results/
```

### Exemplo 2: Gene Típico (~3000 bp → 16384 bp)

```python
# gene.py
gene_seq = open("BRCA1_partial.fasta").read().strip()  # ~3000 bp

# Resize para 16kb
if len(gene_seq) < 16384:
    gene_seq_resized = gene_seq + "N" * (16384 - len(gene_seq))

with open("BRCA1_16kb.fasta", "w") as f:
    f.write(">BRCA1_padded_16kb\n")
    f.write(gene_seq_resized)
```

```bash
python gene.py
python neural_module.py -i BRCA1_16kb.fasta -k API_KEY -o results/
```

### Exemplo 3: Região Grande (~100kb → 131kb)

```bash
# Extrair região de 131kb do genoma
samtools faidx genome.fa chr1:1000000-1131072 > region_131kb.fasta

# Analisar
python neural_module.py -i region_131kb.fasta -k API_KEY -o results/
```

## 🔧 Script Helper para Resize Automático

Criei um script helper que faz o resize automaticamente:

```bash
# resize_fasta.py (criar este script)
```

Uso:
```bash
python resize_fasta.py input.fasta output.fasta --target 2048
```

## 💡 Recomendações

### Para Análise de Variantes
- Use **2048 bp** ou **16384 bp**
- Centre a variante na sequência

### Para Análise de Genes
- Gene pequeno (<10kb): use **16384 bp**
- Gene grande: use **131072 bp** ou **524288 bp**

### Para Análise de Regiões Regulatórias
- Promotor + enhancers próximos: **16384 bp**
- Múltiplos enhancers: **131072 bp**

### Para Análise de Domínios TAD
- Use **524288 bp** ou **1048576 bp**
- Ideal para CONTACT_MAPS

## ⚠️ Notas Importantes

1. **Padding com N**: As bases N não afetarão as predições das bases reais
2. **Centralizar**: Ao fazer padding, considere centralizar sua sequência de interesse
3. **Performance**: Sequências maiores levam mais tempo para processar
4. **Custo**: Sequências maiores podem consumir mais da sua quota de API

## 🆘 Troubleshooting

### Erro: "Sequence length X not supported"
**Solução**: Use um dos 5 tamanhos suportados listados acima

### Sequência exatamente 2048 bp mas ainda dá erro
**Possível causa**: Caracteres invisíveis (espaços, line breaks)
**Solução**: 
```python
seq = seq.replace('\n', '').replace(' ', '').strip()
```

### Como saber o tamanho da minha sequência?
```bash
# Via linha de comando
grep -v ">" sequence.fasta | tr -d '\n' | wc -c

# Via Python
from Bio import SeqIO
rec = SeqIO.read("sequence.fasta", "fasta")
print(len(rec.seq))
```

## 📚 Recursos

- **Verificar tamanhos**: `python check_alphagenome_outputs.py`
- **Documentação**: `OUTPUTS_DISPONIVEIS.md`
- **Exemplo válido**: `example_sequence.fasta` (2048 bp)

---

**Atualizado**: O `neural_module.py` agora valida automaticamente os tamanhos e fornece mensagens claras sobre quais tamanhos são suportados.

*Última atualização: Outubro 2025*

