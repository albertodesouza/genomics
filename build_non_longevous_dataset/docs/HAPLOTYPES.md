# Understanding Haplotypes and Consensus Sequences

## What are Haplotypes?

A **haplotype** is each of the two copies of a chromosome that you inherit - one from each parent (mother and father).

### Why are they important?

All humans are **diploid**: we have 2 copies of each chromosome (except X/Y in males). This means that for any position in the genome, you have:
- **H1 (Haplotype 1)**: The sequence you inherited from one parent
- **H2 (Haplotype 2)**: The sequence you inherited from the other parent

## Consensus Sequences

A **consensus sequence** is generated when you:
1. Take the reference sequence (GRCh38 genome)
2. Apply all individual-specific variants
3. Generate the "real" sequence that individual has

## Practical Example

Imagine a region on chromosome 2:

```
Position:       12345    12350    12355
                  ↓        ↓        ↓
Reference:       A        G        C
Individual HG00096:
  - Variant 1: position 12345, H1=A, H2=T  (heterozygous)
  - Variant 2: position 12350, H1=T, H2=T  (homozygous alternate)
  - Variant 3: position 12355, H1=C, H2=C  (homozygous reference)

Consensus H1:    A T C  (uses haplotype 1 variants)
Consensus H2:    T T C  (uses haplotype 2 variants)
```

## In the Pipeline Context

When `build_window_and_predict.py` runs:

1. **Extracts the reference region** (1 Mb around gene or SNP)
2. **Finds all variants** for that individual in that region from the VCF
3. **Applies variants separately for each haplotype**:
   - `HG00096.H1.window.fixed.fa` ← sequence with haplotype 1 variants
   - `HG00096.H2.window.fixed.fa` ← sequence with haplotype 2 variants

## Why Process Both Haplotypes?

1. **Allele-specific expression**: One haplotype may be expressed differently than the other
2. **Cis effects**: Regulatory variants may affect only one allele
3. **More accurate predictions**: AlphaGenome can predict different activity for H1 vs H2
4. **Imprinting analysis**: Some genes are expressed from only one parental haplotype

## Phasing

1000 Genomes VCFs are **phased**, meaning they know which variant came from which haplotype. This is represented in VCF as:

```
# Unphased (we don't know which allele is on which haplotype)
GT: 0/1    (heterozygous, but phase unknown)

# Phased (we know the phase)
GT: 0|1    (allele 0 on H1, allele 1 on H2)
GT: 1|0    (allele 1 on H1, allele 0 on H2)
```

The pipe symbol `|` indicates phased genotypes, while slash `/` indicates unphased.

## Process Visualization

```
1000 Genomes VCF (phased)
         ↓
    bcftools consensus -H 1    →  HG00096.H1.window.raw.fa
    bcftools consensus -H 2    →  HG00096.H2.window.raw.fa
         ↓
    Adjust to exactly 1 Mb
         ↓
    HG00096.H1.window.fixed.fa  (Final FASTA, 1,000,000 bp)
    HG00096.H2.window.fixed.fa  (Final FASTA, 1,000,000 bp)
         ↓
    AlphaGenome predictions
         ↓
    predictions_H1/  (predictions for haplotype 1)
    predictions_H2/  (predictions for haplotype 2)
```

## Differences You May Encounter

Depending on the variants, H1 and H2 may have:
- Different chromatin accessibility levels (ATAC-seq)
- Different histone modification patterns
- Different predicted gene expression levels
- Different transcription factor binding sites

This is especially relevant for **ancestry analysis**, where one haplotype may come from one ancestry and the other from a different ancestry!

## Technical Details

### How bcftools consensus Works

The `bcftools consensus` command with `-H` flag specifies which haplotype to use:

```bash
# Generate haplotype 1 consensus
bcftools consensus -H 1 -f reference.fa variants.vcf.gz > h1.fa

# Generate haplotype 2 consensus
bcftools consensus -H 2 -f reference.fa variants.vcf.gz > h2.fa
```

For each variant in the VCF:
- If genotype is `0|0`: both haplotypes keep reference allele
- If genotype is `0|1`: H1 keeps reference, H2 gets alternate
- If genotype is `1|0`: H1 gets alternate, H2 keeps reference
- If genotype is `1|1`: both haplotypes get alternate allele

### Handling Complex Variants

- **SNVs (Single Nucleotide Variants)**: Simple substitution
- **Insertions**: Sequence is added to the consensus
- **Deletions**: Sequence is removed from the consensus
- **MNPs (Multiple Nucleotide Polymorphisms)**: Multiple consecutive substitutions

## Ancestry and Haplotypes

In ancestry analysis, haplotypes are particularly interesting because:

1. **Parental origin**: Each haplotype traces back to a specific ancestral lineage
2. **Admixture events**: Recent admixture can create haplotypes with different ancestries
3. **Recombination**: Within a chromosome, segments may have different ancestral origins due to historical recombination

### Example in AISNP Analysis

For an individual with recent European-African admixture:

```
Chromosome 2 (AISNP rs10497191):
  H1: European ancestry → European-specific variants
  H2: African ancestry → African-specific variants

AlphaGenome predictions may show:
  H1: Regulatory pattern typical of European populations
  H2: Regulatory pattern typical of African populations
```

This allows you to:
- Identify ancestry-specific regulatory differences
- Study how genetic background affects gene regulation
- Understand population-specific disease risks

## See Also

- [Main README](../README.md)
- [AISNP Mode Documentation](AISNP_MODE.md)
- [AlphaGenome Predictions Guide](ALPHAGENOME_PREDICTIONS.md)

## Further Reading

- **1000 Genomes phasing methods**: [Nature paper](https://www.nature.com/articles/nature15393)
- **Allele-specific expression**: How variants affect gene expression from each haplotype
- **Haplotype inference**: Statistical methods for determining phase from unphased data
- **Long-read sequencing**: Technologies that can directly sequence entire haplotypes

