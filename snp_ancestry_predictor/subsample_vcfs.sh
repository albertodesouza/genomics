#!/bin/bash
source scripts/start_genomics_universal.sh
mkdir -p /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_subsampled

for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX; do
    input="/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.${chrom}.filtered.SNV_INDEL_SV_phased_panel"
    if [ "$chrom" = "chrX" ]; then
        input="${input}.v2"
    fi
    input="${input}.vcf.gz"
    
    output="/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_subsampled/${chrom}.vcf.gz"
    
    if [ -f "$input" ]; then
        echo "Subsampling $chrom..."
        bcftools view -R random_regions.bed -O z -o "$output" "$input"
        bcftools index -t "$output"
    else
        echo "File $input not found"
    fi
done
echo "Done!"
