# (opcional) Extrair a lista de amostras direto do YAML, se tiver yq instalado
SAMPLES=$(yq '.samples[].sample_id' config_human_30x_monster.yaml | paste -sd' ' -)

# Caso não use yq, defina manualmente:
# SAMPLES="NA12878 NA12891 NA12892"

CUR_DIR=`pwd`
cd /dados/GENOMICS_DATA/top3
echo "current dir: `pwd`"

REF=refs/reference.fa
mkdir -p fasta

# Garanta que a referência e VCFs estão indexados
echo "samtools faidx ${REF}"
samtools faidx ${REF}

for s in $SAMPLES; do
  VCF="vcf/${s}.vcf.gz"
  echo "VCF file: $VCF"
  
  echo "bcftools index -ft ${VCF}"
  bcftools index -ft ${VCF}

  # Dica: seu YAML já configura o bcftools para dividir multialélicos (split). :contentReference[oaicite:6]{index=6}
  # Se precisar normalizar, rode antes:
  # bcftools norm -f "$REF" -m -both -Oz -o "vcf/${s}.norm.vcf.gz" "$VCF" && bcftools index -ft "vcf/${s}.norm.vcf.gz" && VCF="vcf/${s}.norm.vcf.gz"

  # Consenso diploide (respeita GT; se houver phasing, pode escolher haplótipos - ver seção 3)
  echo "bcftools consensus -f ${REF} ${VCF} > fasta/${s}.consensus.fa"
  bcftools consensus -f ${REF} ${VCF} > fasta/${s}.consensus.fa
done

cd $CUR_DIR
