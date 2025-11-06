#!/bin/bash
################################################################################
# Extract AISNPs from 1000 Genomes Project Data
################################################################################
#
# This script downloads 1000 Genomes Project Phase 3 data and extracts
# the 55 Ancestry Informative SNPs (AISNPs) for use with FROGAncestryCalc.
#
# Usage:
#   ./extract_snps_from_1000genomes.sh [options]
#
# Options:
#   -d DIR    Download directory (default: ./1000genomes_data)
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
#   # Extract all samples
#   ./extract_snps_from_1000genomes.sh
#
#   # Extract specific samples
#   echo -e "HG02561\nHG02562\nNA18501" > my_samples.txt
#   ./extract_snps_from_1000genomes.sh -s my_samples.txt -o input/my_data.txt
#
################################################################################

set -e  # Exit on error

# Default parameters
DOWNLOAD_DIR="./1000genomes_data"
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

# 1000 Genomes FTP base URL
FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

echo -e "${YELLOW}Step 1: Downloading 1000 Genomes VCF files...${NC}"
echo "This may take a while (several GB of data)"

# Download VCFs for all chromosomes
for chr in {1..22} X; do
    VCF_FILE="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    TBI_FILE="${VCF_FILE}.tbi"
    
    if [ -f "$VCF_FILE" ]; then
        echo -e "${GREEN}✓ chr${chr} already downloaded${NC}"
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

echo -e "\n${YELLOW}Step 2: Extracting 55 AISNPs from VCF files...${NC}"

# Concatenate and extract SNPs
bcftools concat ALL.chr*.vcf.gz | \
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
    rm -rf "$DOWNLOAD_DIR"/ALL.chr*.vcf.gz*
    echo -e "${GREEN}✓ Cleaned up${NC}"
else
    echo -e "\n${BLUE}ℹ Keeping downloaded VCF files in ${DOWNLOAD_DIR}${NC}"
fi

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}✅ Extraction complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Output file: ${OUTPUT_FILE}"
echo -e "Ready for analysis with: ./run.sh\n"

