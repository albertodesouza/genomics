#!/bin/bash
################################################################################
# Extract AISNPs from Whole Genome Sequencing Data
################################################################################
#
# This script processes WGS data (FASTQ, BAM, or VCF) and extracts the
# 55 Ancestry Informative SNPs for use with FROGAncestryCalc.
#
# Usage:
#   ./extract_snps_from_wgs.sh -i INPUT -t TYPE -o OUTPUT [options]
#
# Required Arguments:
#   -i FILE   Input file (FASTQ, BAM, or VCF)
#   -t TYPE   Input type: fastq, bam, or vcf
#   -o FILE   Output file for FROGAncestryCalc
#
# Optional Arguments:
#   -r FILE   Reference genome (required for fastq/bam)
#   -n NAME   Sample name (default: derived from input file)
#   -2 FILE   Second FASTQ file for paired-end (only for fastq type)
#   -k        Keep intermediate files
#   -h        Show this help message
#
# Requirements:
#   - bcftools, samtools (for all types)
#   - bwa (for fastq type)
#   - Python 3
#
# Examples:
#   # From VCF
#   ./extract_snps_from_wgs.sh -i sample.vcf.gz -t vcf -o input/sample.txt
#
#   # From BAM
#   ./extract_snps_from_wgs.sh -i sample.bam -t bam -r GRCh38.fa -o input/sample.txt
#
#   # From FASTQ (paired-end)
#   ./extract_snps_from_wgs.sh -i sample_R1.fastq.gz -2 sample_R2.fastq.gz \
#       -t fastq -r GRCh38.fa -o input/sample.txt
#
################################################################################

set -e

# Default parameters
INPUT_FILE=""
INPUT_TYPE=""
OUTPUT_FILE=""
REFERENCE=""
SAMPLE_NAME=""
FASTQ2=""
KEEP_FILES=false
SNP_LIST="tools/aisnps_55_list.txt"
THREADS=8

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

show_help() {
    head -n 40 "$0" | grep "^#" | sed 's/^# \?//'
    exit 0
}

# Parse arguments
while getopts "i:t:o:r:n:2:kh" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        t) INPUT_TYPE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        n) SAMPLE_NAME="$OPTARG" ;;
        2) FASTQ2="$OPTARG" ;;
        k) KEEP_FILES=true ;;
        h) show_help ;;
        \?) echo "Invalid option -$OPTARG" >&2; show_help ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_FILE" ] || [ -z "$INPUT_TYPE" ] || [ -z "$OUTPUT_FILE" ]; then
    echo -e "${RED}Error: Missing required arguments${NC}"
    show_help
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

if [[ ! "$INPUT_TYPE" =~ ^(fastq|bam|vcf)$ ]]; then
    echo -e "${RED}Error: Invalid input type. Must be: fastq, bam, or vcf${NC}"
    exit 1
fi

if [[ "$INPUT_TYPE" =~ ^(fastq|bam)$ ]] && [ -z "$REFERENCE" ]; then
    echo -e "${RED}Error: Reference genome required for ${INPUT_TYPE} input${NC}"
    exit 1
fi

# Derive sample name if not provided
if [ -z "$SAMPLE_NAME" ]; then
    SAMPLE_NAME=$(basename "$INPUT_FILE" | sed 's/\.[^.]*$//' | sed 's/\.[^.]*$//')
fi

WORK_DIR="tmp_${SAMPLE_NAME}"
mkdir -p "$WORK_DIR"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}WGS AISNP Extraction Pipeline${NC}"
echo -e "${BLUE}========================================${NC}\n"
echo -e "Sample: ${SAMPLE_NAME}"
echo -e "Input: ${INPUT_FILE}"
echo -e "Type: ${INPUT_TYPE}"
echo -e "Output: ${OUTPUT_FILE}\n"

# Check tools
echo -e "${YELLOW}Checking requirements...${NC}"

check_tool() {
    if ! command -v $1 &> /dev/null; then
        echo -e "${RED}Error: $1 not found${NC}"
        echo "Install with: conda install -c bioconda $1"
        exit 1
    fi
}

check_tool bcftools
check_tool samtools
check_tool python3

if [ "$INPUT_TYPE" == "fastq" ]; then
    check_tool bwa
fi

echo -e "${GREEN}✓ All tools available${NC}\n"

# Process based on input type
case $INPUT_TYPE in
    vcf)
        echo -e "${YELLOW}Processing VCF file...${NC}"
        FINAL_VCF="$INPUT_FILE"
        ;;
        
    bam)
        echo -e "${YELLOW}Step 1: Calling variants from BAM...${NC}"
        bcftools mpileup -f "$REFERENCE" "$INPUT_FILE" | \
            bcftools call -mv -Oz -o "${WORK_DIR}/${SAMPLE_NAME}.vcf.gz"
        bcftools index "${WORK_DIR}/${SAMPLE_NAME}.vcf.gz"
        FINAL_VCF="${WORK_DIR}/${SAMPLE_NAME}.vcf.gz"
        echo -e "${GREEN}✓ Variant calling complete${NC}"
        ;;
        
    fastq)
        echo -e "${YELLOW}Step 1: Aligning reads...${NC}"
        if [ -z "$FASTQ2" ]; then
            # Single-end
            bwa mem -t $THREADS "$REFERENCE" "$INPUT_FILE" | \
                samtools sort -@ $THREADS -o "${WORK_DIR}/${SAMPLE_NAME}.bam"
        else
            # Paired-end
            bwa mem -t $THREADS "$REFERENCE" "$INPUT_FILE" "$FASTQ2" | \
                samtools sort -@ $THREADS -o "${WORK_DIR}/${SAMPLE_NAME}.bam"
        fi
        samtools index "${WORK_DIR}/${SAMPLE_NAME}.bam"
        echo -e "${GREEN}✓ Alignment complete${NC}"
        
        echo -e "\n${YELLOW}Step 2: Calling variants...${NC}"
        bcftools mpileup -f "$REFERENCE" "${WORK_DIR}/${SAMPLE_NAME}.bam" | \
            bcftools call -mv -Oz -o "${WORK_DIR}/${SAMPLE_NAME}.vcf.gz"
        bcftools index "${WORK_DIR}/${SAMPLE_NAME}.vcf.gz"
        FINAL_VCF="${WORK_DIR}/${SAMPLE_NAME}.vcf.gz"
        echo -e "${GREEN}✓ Variant calling complete${NC}"
        ;;
esac

# Extract SNPs
echo -e "\n${YELLOW}Extracting 55 AISNPs...${NC}"
bcftools view -i "ID=@${SNP_LIST}" "$FINAL_VCF" -Oz -o "${WORK_DIR}/${SAMPLE_NAME}_55aisnps.vcf.gz"
bcftools index "${WORK_DIR}/${SAMPLE_NAME}_55aisnps.vcf.gz"
echo -e "${GREEN}✓ SNPs extracted${NC}"

# Convert to FROGAncestryCalc format
echo -e "\n${YELLOW}Converting to FROGAncestryCalc format...${NC}"
python3 tools/vcf_to_frog.py \
    "${WORK_DIR}/${SAMPLE_NAME}_55aisnps.vcf.gz" \
    "$SNP_LIST" \
    "$OUTPUT_FILE"

# Clean up
if [ "$KEEP_FILES" = false ]; then
    echo -e "\n${YELLOW}Cleaning up intermediate files...${NC}"
    rm -rf "$WORK_DIR"
    echo -e "${GREEN}✓ Cleaned up${NC}"
else
    echo -e "\n${BLUE}ℹ Keeping intermediate files in ${WORK_DIR}${NC}"
fi

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}✅ Extraction complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Output file: ${OUTPUT_FILE}"
echo -e "Ready for analysis with: ./run.sh\n"

