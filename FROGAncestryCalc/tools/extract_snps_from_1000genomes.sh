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

# set -e  # Exit on error

# Default parameters
DOWNLOAD_DIR="/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes"
OUTPUT_FILE="input/1000genomes_55aisnps.txt"
SAMPLE_FILE=""
SNP_LIST="tools/aisnps_55_list.txt"
SNP_BED="tools/aisnps_55_grch38.bed"

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
while getopts "d:o:s:h" opt; do
    case $opt in
        d) DOWNLOAD_DIR="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        s) SAMPLE_FILE="$OPTARG" ;;
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

# Create BED file if it doesn't exist
if [ ! -f "$SNP_BED" ]; then
    echo "Creating BED file with SNP coordinates..."
    awk 'NR>1 {print "chr"$3"\t"$4-1"\t"$4"\t"$2}' SNPInfo/55_aisnps_alleles.txt > "$SNP_BED"
    echo -e "${GREEN}✓ Created $SNP_BED${NC}"
fi

echo -e "${GREEN}✓ All requirements met${NC}\n"

# Save original directory
ORIGINAL_DIR="$(pwd)"
echo "Original directory: $ORIGINAL_DIR"

# Create download directory
echo "VCF directory: $DOWNLOAD_DIR"
mkdir -p "$DOWNLOAD_DIR" || {
    echo -e "${RED}Error: Could not create directory $DOWNLOAD_DIR${NC}"
    exit 1
}

cd "$DOWNLOAD_DIR" || {
    echo -e "${RED}Error: Could not change to directory $DOWNLOAD_DIR${NC}"
    exit 1
}

echo "Changed to: $(pwd)"

# 1000 Genomes High Coverage FTP base URL
FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

echo -e "${YELLOW}Step 1: Checking for existing VCF files...${NC}"

# Check if VCFs already exist
EXISTING_COUNT=0
for chr in {1..22} X; do
    VCF_FILE="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    VCF_FILE_V2="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    if [ -f "$VCF_FILE" ] || [ -f "$VCF_FILE_V2" ]; then
        EXISTING_COUNT=$((EXISTING_COUNT + 1))
    fi
done

echo "Found $EXISTING_COUNT out of 23 VCF files"

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
        VCF_FILE_V2="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
        TBI_FILE="${VCF_FILE}.tbi"
        
        # Check both possible file names
        if [ -f "$VCF_FILE" ] || [ -f "$VCF_FILE_V2" ]; then
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

# Check if concatenated file already exists (full genome cache)
CONCAT_VCF="all_chromosomes_concatenated.vcf.gz"
if [ -f "$CONCAT_VCF" ]; then
    echo -e "${GREEN}✓ Found cached concatenated VCF (all chromosomes)${NC}"
    
    # Check if index exists, create if missing
    if [ ! -f "${CONCAT_VCF}.csi" ]; then
        echo -e "${YELLOW}⚠ Index missing, creating index...${NC}"
        echo -e "${YELLOW}Processing...${NC}"
        
        (
            while true; do
                echo -n "."
                sleep 5
            done
        ) &
        INDEX_PID=$!
        
        bcftools index "$CONCAT_VCF"
        
        kill $INDEX_PID 2>/dev/null
        wait $INDEX_PID 2>/dev/null
        echo ""
        
        echo -e "${GREEN}✓ Index created${NC}"
    fi
    
    echo -e "${GREEN}✓ Skipping concatenation step${NC}"
else
    # Build list of VCF files to concatenate (handle both naming patterns)
    VCF_LIST=""
    for chr in {1..22} X; do
        VCF_FILE="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        VCF_FILE_V2="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
        
        if [ -f "$VCF_FILE" ]; then
            VCF_LIST="$VCF_LIST $VCF_FILE"
        elif [ -f "$VCF_FILE_V2" ]; then
            VCF_LIST="$VCF_LIST $VCF_FILE_V2"
        else
            echo -e "${RED}Error: VCF for chr${chr} not found${NC}"
            exit 1
        fi
    done

    VCF_COUNT=$(echo $VCF_LIST | wc -w)
    echo "Found $VCF_COUNT VCF files to process"

    # Concatenate chromosomes
    echo -e "\n${BLUE}ℹ Step 2a: Concatenating 23 chromosomes (may take 3-7 minutes)${NC}"
    echo -e "${BLUE}→ Using 16 threads for faster processing${NC}\n"

    echo "Command being executed:"
    echo "  bcftools concat --threads 16 -Oz -o ${CONCAT_VCF} [23 VCF files]"
    echo ""
    echo -e "${YELLOW}Processing... (a dot will appear every 5 seconds)${NC}"

    # Show progress indicator
    (
        while true; do
            echo -n "."
            sleep 5
        done
    ) &
    PROGRESS_PID=$!

    # Run concatenation
    bcftools concat --threads 16 -Oz -o "$CONCAT_VCF" $VCF_LIST

    # Stop progress indicator
    kill $PROGRESS_PID 2>/dev/null
    wait $PROGRESS_PID 2>/dev/null
    echo ""

    if [ ! -f "$CONCAT_VCF" ]; then
        echo -e "${RED}Error: Failed to concatenate VCFs${NC}"
        exit 1
    fi

    echo "Indexing concatenated VCF (may take 1-2 minutes)..."
    echo -e "${YELLOW}Processing...${NC}"

    # Show progress indicator for indexing
    (
        while true; do
            echo -n "."
            sleep 5
        done
    ) &
    INDEX_PID=$!

    bcftools index "$CONCAT_VCF"

    # Stop progress indicator
    kill $INDEX_PID 2>/dev/null
    wait $INDEX_PID 2>/dev/null
    echo ""

    echo -e "${GREEN}✓ Concatenated and cached all chromosomes (~35GB)${NC}"
    echo -e "${BLUE}ℹ This file will be reused in future runs for much faster extraction${NC}"
fi

# Now extract the 55 SNPs from concatenated file
EXTRACTED_VCF="1000genomes_55aisnps.vcf.gz"
if [ -f "$EXTRACTED_VCF" ]; then
    echo -e "${GREEN}✓ Found cached extracted SNPs file${NC}"
    
    # Check if index exists, create if missing
    if [ ! -f "${EXTRACTED_VCF}.csi" ]; then
        echo "Index missing, creating index (fast)..."
        bcftools index "$EXTRACTED_VCF"
        echo -e "${GREEN}✓ Index created${NC}"
    fi
    
    echo -e "${GREEN}✓ Skipping extraction step${NC}"
else
    echo -e "\n${BLUE}→ Step 2b: Extracting 55 SNPs from concatenated file${NC}"
    echo "Command: bcftools view -R ${SNP_BED} -Oz -o ${EXTRACTED_VCF} ${CONCAT_VCF}"
    
    bcftools view -R "${ORIGINAL_DIR}/${SNP_BED}" -Oz -o "$EXTRACTED_VCF" "$CONCAT_VCF"
    
    if [ ! -f "$EXTRACTED_VCF" ]; then
        echo -e "${RED}Error: Failed to extract SNPs${NC}"
        exit 1
    fi
    
    echo "Indexing extracted VCF (fast, ~few seconds)..."
    bcftools index "$EXTRACTED_VCF"
    
    echo -e "${GREEN}✓ Extracted 55 SNPs (~100KB)${NC}"
fi

# Filter samples if requested
if [ -n "$SAMPLE_FILE" ]; then
    echo -e "\n${YELLOW}Step 3: Filtering specific samples...${NC}"
    
    if [ ! -f "${ORIGINAL_DIR}/${SAMPLE_FILE}" ]; then
        echo -e "${RED}Error: Sample file not found: ${SAMPLE_FILE}${NC}"
        exit 1
    fi
    
    SAMPLES=$(cat "${ORIGINAL_DIR}/${SAMPLE_FILE}" | tr '\n' ',' | sed 's/,$//')
    echo "Filtering samples: ${SAMPLES}"
    bcftools view -s "$SAMPLES" "$EXTRACTED_VCF" -Oz -o 1000genomes_55aisnps_filtered.vcf.gz
    bcftools index 1000genomes_55aisnps_filtered.vcf.gz
    
    FINAL_VCF="1000genomes_55aisnps_filtered.vcf.gz"
    echo -e "${GREEN}✓ Filtered samples${NC}"
else
    FINAL_VCF="$EXTRACTED_VCF"
fi

echo -e "\n${YELLOW}Step 4: Converting to FROGAncestryCalc format...${NC}"

cd "$ORIGINAL_DIR"
python3 tools/vcf_to_frog.py \
    "${DOWNLOAD_DIR}/${FINAL_VCF}" \
    "$SNP_LIST" \
    "$OUTPUT_FILE"

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}✅ Extraction complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Output file: ${OUTPUT_FILE}"
echo -e "\n${BLUE}📁 Cached files in ${DOWNLOAD_DIR}:${NC}"
echo -e "   • all_chromosomes_concatenated.vcf.gz (~35GB)"
echo -e "   • 1000genomes_55aisnps.vcf.gz (~100KB)"
echo -e "\n${BLUE}ℹ Next time you run this script, it will be MUCH faster using the cached files!${NC}"
echo -e "Ready for analysis with: ./run.sh\n"

