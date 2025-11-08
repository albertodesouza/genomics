#!/usr/bin/env bash
set -euo pipefail

# ------------------------------
# Configuração
# ------------------------------
ROOT="/dados/GENOMICS_DATA/top3"
#YAML="/home/lume2/genomics/config_human_30x_monster.yaml"
YAML="/home/alberto/genomics/config_human_30x_monster.yaml"
#ROOT="/dados/GENOMICS_DATA/top4"
#YAML="/home/alberto/genomics/config_human_30x_atena.yaml"
REF="${ROOT}/refs/reference.fa"
OUT="${ROOT}/fasta"
TMP="${OUT}/tmp"

# Ferramentas obrigatórias
for tool in bcftools samtools bedtools awk bgzip; do
  command -v "$tool" >/dev/null 2>&1 || { echo "[ERRO] Ferramenta não encontrada: $tool"; exit 1; }
done

# Amostras (do YAML via yq, se existir; senão, define manualmente)
if command -v yq >/dev/null 2>&1; then
  # vira array bash
  read -r -a SAMPLES <<< "$(yq '.samples[].sample_id' "$YAML" | paste -sd' ' -)"
else
  SAMPLES=(NA12878 NA12891 NA12892)
fi

mkdir -p "$OUT" "$TMP"

# ------------------------------
# Pré-processo fixo
# ------------------------------
echo "[*] Índice da referência: $REF"
samtools faidx "$REF"

# BED “nomeado”: precisa da 4ª coluna (nome do gene) para usar -name+
BED_SRC="${ROOT}/genes/genes.bed"
GENE_LIST="${ROOT}/genes/gene_list.txt"
NAMED_BED="$TMP/genes.named.bed"

if [[ ! -f "$BED_SRC" ]]; then
  echo "[ERRO] BED de genes não encontrado: $BED_SRC"; exit 1
fi

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

echo "[*] BED nomeado para genes: $NAMED_BED"

# ------------------------------
# Loop por amostra
# ------------------------------
for s in "${SAMPLES[@]}"; do
  echo "==> Amostra: $s"
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

  echo "   [OK] Saídas:"
  echo "        - $(basename "$CONS_OUT")"
  echo "        - $(basename "$GENE_OUT")"
done

echo "[✓] Pronto! Todas as saídas estão em: $OUT"
echo "[i] Intermediários ficaram em: $TMP (pode apagar quando quiser)"
