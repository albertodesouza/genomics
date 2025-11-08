#!/usr/bin/env bash
set -euo pipefail

# ========================================
# Generate Gene FASTA Files for Trio Analysis
# ========================================
# This script generates consensus FASTA files per gene from VCF files.
# It processes father, mother, and child samples to create the input
# files needed by genes_difference_count.
#
# Output format: >SampleID|GeneID::chr:start-end(strand)
# ========================================

# ------------------------------
# Configuration - EDIT THESE PATHS
# ------------------------------
# Root directory containing your genomics data
ROOT="/dados/GENOMICS_DATA/top3"

# Path to your configuration YAML file
# Uncomment and modify as needed:
#YAML="/home/lume2/genomics/configs/config_human_30x_monster.yaml"
YAML="/home/alberto/genomics/configs/config_human_30x_monster.yaml"
#ROOT="/dados/GENOMICS_DATA/top4"
#YAML="/home/alberto/genomics/configs/config_human_30x_atena.yaml"

# Reference genome FASTA
REF="${ROOT}/refs/reference.fa"

# Output directory for generated FASTA files
OUT="${ROOT}/fasta"

# Temporary directory for intermediate files
TMP="${OUT}/tmp"

# ------------------------------
# Required Tools Check
# ------------------------------
echo "🔍 Checking required tools..."
for tool in bcftools samtools bedtools awk bgzip; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] Required tool not found: $tool"
    echo "Install with: conda install -c bioconda $tool"
    exit 1
  fi
done
echo "✓ All required tools are available"

# ------------------------------
# Sample Detection
# ------------------------------
# Try to read samples from YAML (if yq is available), otherwise use defaults
if command -v yq >/dev/null 2>&1 && [[ -f "$YAML" ]]; then
  echo "📄 Reading samples from YAML: $YAML"
  read -r -a SAMPLES <<< "$(yq '.samples[].sample_id' "$YAML" | paste -sd' ' -)"
  echo "   Found samples: ${SAMPLES[*]}"
else
  echo "⚠️  yq not found or YAML missing, using default samples"
  # Default trio samples (father, mother, child)
  SAMPLES=(NA12891 NA12892 NA12878)
  echo "   Using: ${SAMPLES[*]} (edit script to customize)"
fi

mkdir -p "$OUT" "$TMP"

# ------------------------------
# Pre-processing
# ------------------------------
echo ""
echo "📁 Working directories:"
echo "   Input VCFs: ${ROOT}/vcf/"
echo "   Reference: $REF"
echo "   Output: $OUT"
echo "   Temp: $TMP"
echo ""

if [[ ! -f "$REF" ]]; then
  echo "[ERROR] Reference genome not found: $REF"
  exit 1
fi

echo "🔨 Indexing reference genome..."
samtools faidx "$REF"

# ------------------------------
# Gene BED File Preparation
# ------------------------------
# BED file needs 4th column (gene name) for bedtools getfasta -name+
BED_SRC="${ROOT}/genes/genes.bed"
GENE_LIST="${ROOT}/genes/gene_list.txt"
NAMED_BED="$TMP/genes.named.bed"

if [[ ! -f "$BED_SRC" ]]; then
  echo "[ERROR] Gene BED file not found: $BED_SRC"
  echo "Make sure you have run the gene extraction step of the pipeline"
  exit 1
fi

echo "🧬 Preparing gene BED file..."

# Detecta se já tem 4ª coluna; se não tiver, tenta colar gene_list; se não der, usa coords como nome
NF_FIRST_NONCOMMENT=$(awk '($0 !~ /^#/ && NF>0){print NF; exit}' "$BED_SRC" || true)
if [[ "${NF_FIRST_NONCOMMENT:-0}" -ge 4 ]]; then
  awk 'BEGIN{OFS="\t"} $0 ~ /^#/ {next} {print $1,$2,$3,$4}' "$BED_SRC" > "$NAMED_BED"
elif [[ -f "$GENE_LIST" ]] && [[ $(wc -l < "$GENE_LIST") -eq $(grep -vc '^#' "$BED_SRC") ]]; then
  paste "$BED_SRC" "$GENE_LIST" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' > "$NAMED_BED"
else
  # sem gene_list compatível: usa coord como “nome”
  awk 'BEGIN{OFS="\t"} $0 !~ /^#/ {print $1,$2,$3,$1":"$2"-"$3}' "$BED_SRC" > "$NAMED_BED"
fi

echo "✓ Gene BED ready: $NAMED_BED"
echo ""

# ------------------------------
# Process Each Sample
# ------------------------------
echo "🚀 Processing samples..."
echo "========================================"
for s in "${SAMPLES[@]}"; do
  echo ""
  echo "📊 Sample: $s"
  VCF_IN="${ROOT}/vcf/${s}.vcf.gz"
  if [[ ! -f "$VCF_IN" ]]; then
    echo "   [AVISO] VCF ausente: $VCF_IN — pulando."
    continue
  fi

  # Trabalha SEM tocar em vcf/: tudo vai para fasta/tmp/
  VCF_SORTED="$TMP/${s}.sorted.vcf.gz"
  VCF_NORM="$TMP/${s}.norm.vcf.gz"
  VCF_NOOV="$TMP/${s}.norm.noov.vcf.gz"
  VCF_USE=""

  echo "   [*] Ordenando VCF…"
  bcftools sort -Oz -o "$VCF_SORTED" "$VCF_IN"
  bcftools index -ft "$VCF_SORTED"

  echo "   [*] Normalizando (left-align & decompondo multialélicos/MNPs)…"
  bcftools norm -f "$REF" -m -both -Oz -o "$VCF_NORM" "$VCF_SORTED"
  bcftools index -ft "$VCF_NORM"

  echo "   [*] Removendo sobreposições (se plugin disponível)…"
  if bcftools +remove-overlaps "$VCF_NORM" -d -Ou >/dev/null 2>&1; then
    bcftools +remove-overlaps "$VCF_NORM" -d | bgzip -c > "$VCF_NOOV"
    bcftools index -ft "$VCF_NOOV"
    VCF_USE="$VCF_NOOV"
  else
    echo "   [AVISO] Plugin +remove-overlaps indisponível; seguindo com VCF normalizado."
    VCF_USE="$VCF_NORM"
  fi

  # Nome real da amostra dentro do VCF (no seu caso costuma ser bam/<sample>.mkdup.bam)
  SAMPLE_IN_VCF=$(bcftools query -l "$VCF_USE" | head -n1)
  if [[ -z "$SAMPLE_IN_VCF" ]]; then
    echo "   [ERRO] VCF sem coluna de amostra/GT: $VCF_IN — pulando."
    continue
  fi

  # 1) FASTA consenso por cromossomo (gera sem prefixo; etiqueta depois)
  CONS_CORE="$TMP/${s}.consensus.core.fa"
  echo "   [*] Gerando consenso (IUPAC p/ heterozigotos por padrão)…"
  bcftools consensus -f "$REF" -s "$SAMPLE_IN_VCF" "$VCF_USE" > "$CONS_CORE"

  # Etiquetar headers com “Amostra|chr”
  CONS_OUT="$OUT/${s}.consensus.fa"
  awk -v P="${s}|" 'BEGIN{OFS=""} /^>/{sub(/^>/,">" P); print; next} {print}' "$CONS_CORE" > "$CONS_OUT"

  # 2) FASTA por gene (usando o consenso “core” para bater cromossomos com o BED)
  GENE_CORE="$TMP/${s}.genes.core.fa"
  GENE_OUT="$OUT/${s}.genes.consensus.fa"
  echo "   [*] Extraindo sequências por gene (bedtools getfasta -s -name+ )…"
  bedtools getfasta -fi "$CONS_CORE" -bed "$NAMED_BED" -s -name+ -fo "$GENE_CORE"
  # Prefixar com amostra: >Amostra|GENE|chr:start-end(strand)
  awk -v P="${s}|" 'BEGIN{OFS=""} /^>/{sub(/^>/,">" P); print; next} {print}' "$GENE_CORE" > "$GENE_OUT"

  echo "   ✓ Outputs generated:"
  echo "      - $(basename "$CONS_OUT") (chromosome consensus)"
  echo "      - $(basename "$GENE_OUT") (gene consensus - INPUT FOR genes_difference_count)"
done

echo ""
echo "========================================"
echo "✅ COMPLETE! All outputs in: $OUT"
echo ""
echo "📄 Generated files:"
for s in "${SAMPLES[@]}"; do
  if [[ -f "$OUT/${s}.genes.consensus.fa" ]]; then
    size=$(ls -lh "$OUT/${s}.genes.consensus.fa" | awk '{print $5}')
    genes=$(grep -c "^>" "$OUT/${s}.genes.consensus.fa" || echo "0")
    echo "   - ${s}.genes.consensus.fa ($size, $genes genes)"
  fi
done
echo ""
echo "💡 Next steps:"
echo "   1. Edit genes_difference_count.cpp to set these input paths"
echo "   2. Compile: make"
echo "   3. Run: ./genes_difference_count"
echo ""
echo "🗑️  Intermediate files in: $TMP (can be deleted)"
