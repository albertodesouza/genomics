#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------
# Configuração
# ------------------------------------
ROOT="/dados/GENOMICS_DATA/top3"
COMP="${ROOT}/comparisons"
OUTDIR="${COMP}"
OUTCSV="${OUTDIR}/pairwise_diff_summary.csv"

# Checagens
for t in bcftools awk; do
  command -v "$t" >/dev/null 2>&1 || { echo "[ERRO] Ferramenta não encontrada: $t"; exit 1; }
done
[[ -d "$COMP" ]] || { echo "[ERRO] Pasta não encontrada: $COMP"; exit 1; }

echo "pair,mode,total_records,both_called,matches,mismatches,missing_one,both_missing,mismatch_rate" > "$OUTCSV"

pair_stats () {
  local VCF="$1"
  local MODE="$2"   # "all" ou "snv"

  # Amostras do VCF (esperado 2 colunas)
  local PAIR
  PAIR="$(bcftools query -l "$VCF" | paste -sd'_' -)"

  # Stream de linhas: CHROM POS GT1 GT2 (apenas PASS)
  if [[ "$MODE" == "snv" ]]; then
    # SNVs bialélicos com PASS
    STREAM_CMD="bcftools view -v snps -m2 -M2 -i 'FILTER==\"PASS\"' \"$VCF\" | bcftools query -f '%CHROM\t%POS[\t%GT]\n' -"
  else
    # Todas variantes com PASS
    STREAM_CMD="bcftools query -i 'FILTER==\"PASS\"' -f '%CHROM\t%POS[\t%GT]\n' \"$VCF\""
  fi

  # Conta
  IFS=, read -r TOT BOTH MATCH MISM MISS1 RATE BOTHMISS < <(
    bash -c "$STREAM_CMD" | \
    awk 'BEGIN{FS="\t"}
         {
           gt1=$3; gt2=$4;
           tot++;
           gsub(/\|/,"/",gt1); gsub(/\|/,"/",gt2);          # ignora fase
           miss1 = (gt1 ~ /\./); miss2 = (gt2 ~ /\./);

           if (miss1 && miss2) { both_missing++ }
           else if (miss1 || miss2) { missing_one++ }
           else {
             # normaliza diplóide para ordem canônica (0/1 == 1/0)
             n1=split(gt1,a,"/"); n2=split(gt2,b,"/");
             if (n1==2 && a[1]+0 > a[2]+0) { t=a[1]; a[1]=a[2]; a[2]=t }
             if (n2==2 && b[1]+0 > b[2]+0) { t=b[1]; b[1]=b[2]; b[2]=t }
             if (n1==2) gt1=a[1] "/" a[2];
             if (n2==2) gt2=b[1] "/" b[2];

             both_called++;
             if (gt1==gt2) matches++;
             else mismatches++;
           }
         }
         END{
           rate = (both_called>0) ? mismatches / both_called : 0;
           printf("%d,%d,%d,%d,%d,%.6f,%d\n",
                 tot, both_called, matches, mismatches, missing_one, rate, both_missing);
         }'
  )

  echo "${PAIR},${MODE},${TOT},${BOTH},${MATCH},${MISM},${MISS1},${BOTHMISS},${RATE}" >> "$OUTCSV"

  # Saída amigável no terminal
  printf "\n[%s] mode=%s\n" "$PAIR" "$MODE"
  printf "  total records (PASS):     %s\n" "$TOT"
  printf "  both called:              %s\n" "$BOTH"
  printf "    ├─ matches (same GT):   %s\n" "$MATCH"
  printf "    └─ mismatches:          %s\n" "$MISM"
  printf "  missing (one sample):     %s\n" "$MISS1"
  printf "  missing (both samples):   %s\n" "$BOTHMISS"
  printf "  mismatch rate (both-called): %.6f\n" "$RATE"
}

# Loop sobre os VCFs pareados
shopt -s nullglob
PAIR_VCFS=( "${COMP}"/*_vs_*.merge.vcf.gz )
if [[ ${#PAIR_VCFS[@]} -eq 0 ]]; then
  echo "[ERRO] Nenhum arquivo *_vs_*.merge.vcf.gz encontrado em ${COMP}"
  exit 2
fi

for VCF in "${PAIR_VCFS[@]}"; do
  echo "[*] Processando: $(basename "$VCF")"
  pair_stats "$VCF" "all"
  pair_stats "$VCF" "snv"
done

echo -e "\n[OK] Resumo salvo em: ${OUTCSV}"
