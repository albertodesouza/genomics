#!/bin/bash
################################################################################
# Extract AISNPs from 1000 Genomes Project Data
################################################################################
#
# This script downloads 1000 Genomes High Coverage (GRCh38) data and extracts
# the 55 Ancestry Informative SNPs (AISNPs) for use with FROGAncestryCalc.
# If VCFs already exist in the specified directory, download is skipped.
#
# Usage:
#   ./extract_snps_from_1000genomes.sh [options]
#
# Options:
#   -d DIR    VCF directory (default: /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes)
#             Script will auto-detect existing VCFs and skip download if present
#   -o FILE   Output file (default: input/1000genomes_55aisnps.txt)
#   -s FILE   Sample list (optional, extracts specific samples only)
#   -k        Keep downloaded VCF files (default: delete after extraction)
#   -h        Show this help message
#
# Requirements:
#   - bcftools (install via: conda install -c bioconda bcftools)
#   - wget
#   - Python 3
#
# Example:
#   # Extract all samples (uses existing VCFs if available)
#   ./extract_snps_from_1000genomes.sh
#
#   # Extract specific samples
#   echo -e "HG02561\nHG02562\nHG03055" > my_samples.txt
#   ./extract_snps_from_1000genomes.sh -s my_samples.txt -o input/my_data.txt
#
################################################################################

set -e  # Exit on error

# Default parameters
DOWNLOAD_DIR="/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes"
OUTPUT_FILE="input/1000genomes_55aisnps.txt"
SAMPLE_FILE=""
KEEP_FILES=false
SNP_LIST="tools/aisnps_55_list.txt"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Help message
show_help() {
    head -n 30 "$0" | grep "^#" | sed 's/^# \?//'
    exit 0
}

# Parse command line arguments
while getopts "d:o:s:kh" opt; do
    case $opt in
        d) DOWNLOAD_DIR="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        s) SAMPLE_FILE="$OPTARG" ;;
        k) KEEP_FILES=true ;;
        h) show_help ;;
        \?) echo "Invalid option -$OPTARG" >&2; show_help ;;
    esac
done

# Print header
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}1000 Genomes AISNP Extraction Pipeline${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Check requirements
echo -e "${YELLOW}Checking requirements...${NC}"

if ! command -v bcftools &> /dev/null; then
    echo -e "${RED}Error: bcftools not found${NC}"
    echo "Install with: conda install -c bioconda bcftools"
    exit 1
fi

if ! command -v wget &> /dev/null; then
    echo -e "${RED}Error: wget not found${NC}"
    echo "Install with: sudo apt-get install wget"
    exit 1
fi

if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: python3 not found${NC}"
    exit 1
fi

if [ ! -f "$SNP_LIST" ]; then
    echo -e "${RED}Error: SNP list not found: $SNP_LIST${NC}"
    echo "Make sure you're running this from the FROGAncestryCalc directory"
    exit 1
fi

echo -e "${GREEN}✓ All requirements met${NC}\n"

# Create download directory
mkdir -p "$DOWNLOAD_DIR"
cd "$DOWNLOAD_DIR"

# 1000 Genomes High Coverage FTP base URL
FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

echo -e "${YELLOW}Step 1: Checking for existing VCF files...${NC}"

# Check if VCFs already exist
EXISTING_COUNT=0
for chr in {1..22} X; do
    VCF_FILE="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    if [ -f "$VCF_FILE" ]; then
        ((EXISTING_COUNT++))
    fi
done

if [ $EXISTING_COUNT -eq 23 ]; then
    echo -e "${GREEN}✓ All 23 VCF files found in ${DOWNLOAD_DIR}${NC}"
    echo -e "${GREEN}✓ Skipping download step${NC}"
elif [ $EXISTING_COUNT -gt 0 ]; then
    echo -e "${YELLOW}⚠ Found ${EXISTING_COUNT}/23 VCF files${NC}"
    echo -e "${YELLOW}→ Downloading missing chromosomes...${NC}"
    NEED_DOWNLOAD=true
else
    echo -e "${YELLOW}→ No existing VCFs found${NC}"
    echo -e "${YELLOW}→ Downloading 1000 Genomes High Coverage VCF files...${NC}"
    echo "This may take a while (several GB of data)"
    NEED_DOWNLOAD=true
fi

# Download VCFs if needed
if [ "$NEED_DOWNLOAD" = true ]; then
    for chr in {1..22} X; do
        VCF_FILE="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        TBI_FILE="${VCF_FILE}.tbi"
        
        if [ -f "$VCF_FILE" ]; then
            echo -e "${GREEN}✓ chr${chr} already present${NC}"
        else
            echo -e "${BLUE}→ Downloading chromosome ${chr}...${NC}"
            wget -q --show-progress "${FTP_BASE}/${VCF_FILE}" || {
                echo -e "${RED}Failed to download chr${chr}${NC}"
                exit 1
            }
            wget -q "${FTP_BASE}/${TBI_FILE}" || {
                echo -e "${RED}Failed to download chr${chr} index${NC}"
                exit 1
            }
        fi
    done
fi

echo -e "\n${YELLOW}Step 2: Extracting 55 AISNPs from VCF files...${NC}"

# Concatenate and extract SNPs
bcftools concat 1kGP_high_coverage_Illumina.chr*.vcf.gz | \
    bcftools view -i "ID=@../${SNP_LIST}" -Oz -o 1000genomes_55aisnps.vcf.gz

bcftools index 1000genomes_55aisnps.vcf.gz

echo -e "${GREEN}✓ Extracted SNPs${NC}"

# Filter samples if requested
if [ -n "$SAMPLE_FILE" ]; then
    echo -e "\n${YELLOW}Step 3: Filtering specific samples...${NC}"
    
    if [ ! -f "../${SAMPLE_FILE}" ]; then
        echo -e "${RED}Error: Sample file not found: ${SAMPLE_FILE}${NC}"
        exit 1
    fi
    
    SAMPLES=$(cat "../${SAMPLE_FILE}" | tr '\n' ',' | sed 's/,$//')
    bcftools view -s "$SAMPLES" 1000genomes_55aisnps.vcf.gz -Oz -o 1000genomes_55aisnps_filtered.vcf.gz
    bcftools index 1000genomes_55aisnps_filtered.vcf.gz
    
    FINAL_VCF="1000genomes_55aisnps_filtered.vcf.gz"
else
    FINAL_VCF="1000genomes_55aisnps.vcf.gz"
fi

echo -e "\n${YELLOW}Step 4: Converting to FROGAncestryCalc format...${NC}"

cd ..
python3 tools/vcf_to_frog.py \
    "${DOWNLOAD_DIR}/${FINAL_VCF}" \
    "$SNP_LIST" \
    "$OUTPUT_FILE"

# Clean up if requested
if [ "$KEEP_FILES" = false ]; then
    echo -e "\n${YELLOW}Cleaning up downloaded files...${NC}"
    rm -rf "$DOWNLOAD_DIR"/1kGP_high_coverage_Illumina.chr*.vcf.gz*
    echo -e "${GREEN}✓ Cleaned up${NC}"
else
    echo -e "\n${BLUE}ℹ Keeping VCF files in ${DOWNLOAD_DIR}${NC}"
fi

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}✅ Extraction complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Output file: ${OUTPUT_FILE}"
echo -e "Ready for analysis with: ./run.sh\n"

