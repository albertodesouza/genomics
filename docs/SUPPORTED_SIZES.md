# ⚠️ Sequence Sizes Supported by AlphaGenome

## 🎯 Important Information

AlphaGenome **only accepts sequences with specific sizes**:

| Size | Bytes | Recommended Use |
|------|-------|-----------------|
| **2,048 bp** | 2 KB | Small regions, quick tests |
| **16,384 bp** | 16 KB | Individual genes |
| **131,072 bp** | 128 KB | Multiple genes, clusters |
| **524,288 bp** | 512 KB | Large genomic regions |
| **1,048,576 bp** | 1 MB | Extensive chromosomal regions |

## ❌ What Happens with Other Sizes?

If you try to use a sequence of 600 bp, 5000 bp, or any other unlisted size, you'll receive the error:

```
ValueError: Sequence length 600 not supported by the model. 
Supported lengths: [2048, 16384, 131072, 524288, 1048576]
```

## ✅ How to Fix

### Option 1: Pad (Add Bases)

If your sequence is **smaller** than the minimum size (2048 bp), add N bases:

```python
sequence = "ATCGATCG..."  # 600 bp
padded = sequence + "N" * (2048 - len(sequence))  # Now has 2048 bp
```

### Option 2: Truncate

If your sequence is larger but not a supported size, truncate to the next smaller size:

```python
sequence = "ATCGATCG..." * 1000  # 5000 bp
truncated = sequence[:2048]  # Use first 2048 bp
```

### Option 3: Resize to Next Supported Size

```python
def resize_to_supported(seq):
    supported = [2048, 16384, 131072, 524288, 1048576]
    
    # Find the nearest supported size
    for size in supported:
        if len(seq) <= size:
            if len(seq) < size:
                # Add Ns
                return seq + "N" * (size - len(seq))
            return seq
    
    # If larger than 1MB, truncate
    return seq[:1048576]
```

## 📝 Practical Examples

### Example 1: Small Sequence (600 bp → 2048 bp)

```bash
# Create sequence with padding in Python
python3 << 'EOF'
seq = "ATCG" * 150  # 600 bp
padded = seq + "N" * (2048 - len(seq))
with open("seq_2kb.fasta", "w") as f:
    f.write(">my_sequence_2kb\n")
    f.write(padded + "\n")
EOF

# Analyze
python neural_module/neural_module.py -i seq_2kb.fasta -k API_KEY -o results/
```

### Example 2: Typical Gene (~3000 bp → 16384 bp)

```python
# gene.py
gene_seq = open("BRCA1_partial.fasta").read().strip()  # ~3000 bp

# Resize to 16kb
if len(gene_seq) < 16384:
    gene_seq_resized = gene_seq + "N" * (16384 - len(gene_seq))

with open("BRCA1_16kb.fasta", "w") as f:
    f.write(">BRCA1_padded_16kb\n")
    f.write(gene_seq_resized)
```

```bash
python gene.py
python neural_module/neural_module.py -i BRCA1_16kb.fasta -k API_KEY -o results/
```

### Example 3: Large Region (~100kb → 131kb)

```bash
# Extract 131kb region from genome
samtools faidx genome.fa chr1:1000000-1131072 > region_131kb.fasta

# Analyze
python neural_module/neural_module.py -i region_131kb.fasta -k API_KEY -o results/
```

## 🔧 Helper Script for Automatic Resize

Created a helper script that does the resize automatically:

```bash
# resize_fasta.py (create this script)
```

Usage:
```bash
python resize_fasta.py input.fasta output.fasta --target 2048
```

## 💡 Recommendations

### For Variant Analysis
- Use **2048 bp** or **16384 bp**
- Center the variant in the sequence

### For Gene Analysis
- Small gene (<10kb): use **16384 bp**
- Large gene: use **131072 bp** or **524288 bp**

### For Regulatory Region Analysis
- Promoter + nearby enhancers: **16384 bp**
- Multiple enhancers: **131072 bp**

### For TAD Domain Analysis
- Use **524288 bp** or **1048576 bp**
- Ideal for CONTACT_MAPS

## ⚠️ Important Notes

1. **Padding with N**: N bases won't affect predictions of real bases
2. **Centering**: When padding, consider centering your sequence of interest
3. **Performance**: Larger sequences take longer to process
4. **Cost**: Larger sequences may consume more of your API quota

## 🆘 Troubleshooting

### Error: "Sequence length X not supported"
**Solution**: Use one of the 5 supported sizes listed above

### Sequence is exactly 2048 bp but still gives error
**Possible cause**: Invisible characters (spaces, line breaks)
**Solution**: 
```python
seq = seq.replace('\n', '').replace(' ', '').strip()
```

### How to check my sequence size?
```bash
# Via command line
grep -v ">" sequence.fasta | tr -d '\n' | wc -c

# Via Python
from Bio import SeqIO
rec = SeqIO.read("sequence.fasta", "fasta")
print(len(rec.seq))
```

## 📚 Resources

- **Check sizes**: `python scripts/check_alphagenome_outputs.py`
- **Documentation**: See `docs/AVAILABLE_OUTPUTS.md`
- **Valid example**: `example_sequence.fasta` (2048 bp)

---

**Updated**: `neural_module.py` now automatically validates sizes and provides clear messages about which sizes are supported.

*Last updated: October 2025*

