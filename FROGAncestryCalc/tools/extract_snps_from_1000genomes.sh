#!/bin/bash
################################################################################
# Extract AISNPs from 1000 Genomes Project Data
################################################################################
#
# This script downloads 1000 Genomes data (Phase 3 GRCh37 or High Coverage GRCh38)
# and extracts the 55 Ancestry Informative SNPs (AISNPs) for use with FROGAncestryCalc.
# If VCFs already exist in the specified directory, download is skipped.
#
# Usage:
#   ./extract_snps_from_1000genomes.sh [options]
#
# Options:
#   -b BUILD  Genome build: grch37 or grch38 (default: grch37)
#   -d DIR    VCF directory (default: auto-set based on build)
#             GRCh37: /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes_GRCh37
#             GRCh38: /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes
#             Script will auto-detect existing VCFs and skip download if present
#   -o FILE   Output file (default: input/1000genomes_55aisnps.txt)
#   -s FILE   Sample list (optional, extracts specific samples only)
#   -h        Show this help message
#
# Requirements:
#   - bcftools (install via: conda install -c bioconda bcftools)
#   - wget
#   - Python 3
#   - For GRCh38: Run tools/convert_grch37_to_grch38.py first to generate coordinates
#
# Examples:
#   # Extract all samples using GRCh37 (default, uses existing VCFs if available)
#   ./extract_snps_from_1000genomes.sh
#
#   # Extract using GRCh38 High Coverage data
#   ./extract_snps_from_1000genomes.sh -b grch38
#
#   # Extract specific samples with GRCh37
#   echo -e "HG02561\nHG02562\nHG03055" > my_samples.txt
#   ./extract_snps_from_1000genomes.sh -b grch37 -s my_samples.txt -o input/my_data.txt
#
################################################################################

# set -e  # Exit on error

# Default parameters
BUILD="GRCh38"
DOWNLOAD_DIR=""  # Will be set based on build
OUTPUT_FILE="input/1000genomes_55aisnps.txt"
SAMPLE_FILE=""
SNP_LIST="tools/aisnps_55_list.txt"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Help message
show_help() {
    head -n 40 "$0" | grep "^#" | sed 's/^# \?//'
    exit 0
}

# Parse command line arguments
while getopts "b:d:o:s:h" opt; do
    case $opt in
        b) BUILD="$OPTARG" ;;
        d) DOWNLOAD_DIR="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        s) SAMPLE_FILE="$OPTARG" ;;
        h) show_help ;;
        \?) echo "Invalid option -$OPTARG" >&2; show_help ;;
    esac
done

# Validate build option
BUILD=$(echo "$BUILD" | tr '[:upper:]' '[:lower:]')
if [[ "$BUILD" != "grch37" && "$BUILD" != "grch38" ]]; then
    echo -e "${RED}Error: Invalid build '$BUILD'. Must be 'grch37' or 'grch38'${NC}"
    exit 1
fi

# Set build-specific defaults
if [ -z "$DOWNLOAD_DIR" ]; then
    if [ "$BUILD" = "grch38" ]; then
        DOWNLOAD_DIR="/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes"
    else
        DOWNLOAD_DIR="/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes_GRCh37"
    fi
fi

# Set build-specific files and URLs
if [ "$BUILD" = "grch38" ]; then
    ALLELES_FILE="SNPInfo/55_aisnps_alleles_grch38.txt"
    SNP_BED="tools/aisnps_55_grch38.bed"
    FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
    VCF_PATTERN="1kGP_high_coverage_Illumina"
    BUILD_DISPLAY="GRCh38/hg38 (High Coverage)"
else
    ALLELES_FILE="SNPInfo/55_aisnps_alleles.txt"
    SNP_BED="tools/aisnps_55_grch37.bed"
    FTP_BASE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
    VCF_PATTERN="phase3"
    BUILD_DISPLAY="GRCh37/hg19 (Phase 3)"
fi

# Print header
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}1000 Genomes AISNP Extraction Pipeline${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Genome Build: ${BUILD_DISPLAY}${NC}"
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

# Check for alleles file
if [ ! -f "$ALLELES_FILE" ]; then
    echo -e "${RED}Error: Alleles file not found: $ALLELES_FILE${NC}"
    if [ "$BUILD" = "grch38" ]; then
        echo -e "${YELLOW}For GRCh38, you need to generate the coordinates file first:${NC}"
        echo "  python3 tools/convert_grch37_to_grch38.py"
    fi
    exit 1
fi

# Create BED file if it doesn't exist
if [ ! -f "$SNP_BED" ]; then
    echo "Creating BED file with SNP coordinates from $ALLELES_FILE..."
    if [ "$BUILD" = "grch37" ]; then
        # GRCh37/Phase 3: SEM prefixo chr (usa apenas 1, 2, 3, ..., X)
        awk 'NR>1 {print $3"\t"$4-1"\t"$4"\t"$2}' "$ALLELES_FILE" > "$SNP_BED"
    else
        # GRCh38: COM prefixo chr (usa chr1, chr2, chr3, ..., chrX)
        awk 'NR>1 {print "chr"$3"\t"$4-1"\t"$4"\t"$2}' "$ALLELES_FILE" > "$SNP_BED"
    fi
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

echo -e "${YELLOW}Step 1: Checking for existing VCF files...${NC}"

# Function to get VCF filename for a chromosome
get_vcf_filename() {
    local chr=$1
    if [ "$BUILD" = "grch38" ]; then
        # High Coverage pattern (might have .v2 suffix for some chromosomes)
        echo "1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel"
    else
        # Phase 3 pattern
        # Note: chrX uses v1c instead of v5b
        if [ "$chr" = "X" ]; then
            echo "ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes"
        else
            echo "ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes"
        fi
    fi
}

# Robust download function with retry and timeout
download_file_robust() {
    local url=$1
    local output=$2
    local max_attempts=5
    local stall_timeout=120  # Kill if no progress for 2 minutes
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        if [ $attempt -eq 1 ]; then
            echo -e "${BLUE}  Starting download...${NC}"
        else
            echo -e "${BLUE}  Retry attempt $attempt/$max_attempts...${NC}"
        fi
        
        # Remove partial file if exists on first attempt
        if [ $attempt -eq 1 ]; then
            [ -f "${output}.partial" ] && rm -f "${output}.partial"
        fi
        
        # Get start time for speed calculation
        start_time=$(date +%s)
        
        # Start download in background with watchdog
        if command -v wget &> /dev/null; then
            # Use wget with progress bar in background
            (
                wget --timeout=60 \
                     --read-timeout=60 \
                     --tries=1 \
                     --continue \
                     --show-progress \
                     --progress=bar:force:noscroll \
                     -O "${output}.partial" \
                     "$url" 2>&1
            ) &
            download_pid=$!
        # Fallback to curl if wget not available
        elif command -v curl &> /dev/null; then
            echo -e "${YELLOW}  Using curl (wget not available)${NC}"
            (
                curl --connect-timeout 30 \
                     --max-time 0 \
                     --speed-time 120 \
                     --speed-limit 1000 \
                     --retry 0 \
                     --progress-bar \
                     --continue-at - \
                     -o "${output}.partial" \
                     "$url"
            ) &
            download_pid=$!
        else
            echo -e "${RED}Error: Neither wget nor curl is available${NC}"
            return 1
        fi
        
        # Watchdog: monitor if file is growing
        last_size=0
        stall_count=0
        max_stalls=3  # Allow 3 checks without progress before killing
        
        while kill -0 $download_pid 2>/dev/null; do
            sleep 20  # Check every 20 seconds
            
            if [ -f "${output}.partial" ]; then
                current_size=$(stat -f%z "${output}.partial" 2>/dev/null || stat -c%s "${output}.partial" 2>/dev/null || echo 0)
                
                if [ "$current_size" -eq "$last_size" ]; then
                    stall_count=$((stall_count + 1))
                    echo -e "${YELLOW}  ⚠ No progress detected (check $stall_count/$max_stalls)${NC}"
                    
                    if [ $stall_count -ge $max_stalls ]; then
                        echo -e "${RED}  ✗ Download stalled - killing process${NC}"
                        kill -9 $download_pid 2>/dev/null
                        wait $download_pid 2>/dev/null
                        break
                    fi
                else
                    # Progress detected, reset counter
                    if [ $stall_count -gt 0 ]; then
                        echo -e "${GREEN}  ✓ Progress resumed${NC}"
                    fi
                    stall_count=0
                    last_size=$current_size
                fi
            fi
        done
        
        # Wait for download to finish
        wait $download_pid 2>/dev/null
        download_status=$?
        
        # Calculate download time
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        
        # Check if download succeeded
        if [ $download_status -eq 0 ] && [ -f "${output}.partial" ]; then
            # Check if file has content (not empty or too small)
            file_size=$(stat -f%z "${output}.partial" 2>/dev/null || stat -c%s "${output}.partial" 2>/dev/null || echo 0)
            
            if [ "$file_size" -gt 1000 ]; then
                # Move to final location
                mv "${output}.partial" "$output"
                
                # Format file size
                if command -v numfmt &> /dev/null; then
                    size_human=$(numfmt --to=iec-i --suffix=B $file_size)
                else
                    size_human=$(echo "scale=2; $file_size / 1024 / 1024 / 1024" | bc 2>/dev/null || echo "$file_size")
                    size_human="${size_human} GB"
                fi
                
                # Calculate average speed
                if [ $duration -gt 0 ]; then
                    speed_bps=$((file_size / duration))
                    if command -v numfmt &> /dev/null; then
                        speed_human=$(numfmt --to=iec-i --suffix=B/s $speed_bps)
                    else
                        speed_mbps=$(echo "scale=2; $speed_bps / 1024 / 1024" | bc 2>/dev/null || echo "?")
                        speed_human="${speed_mbps} MB/s"
                    fi
                    echo -e "${GREEN}  ✓ Download complete: ${size_human} in ${duration}s (avg: ${speed_human})${NC}"
                else
                    echo -e "${GREEN}  ✓ Download complete: ${size_human}${NC}"
                fi
                return 0
            else
                echo -e "${YELLOW}  ⚠ Downloaded file too small ($file_size bytes), retrying...${NC}"
                rm -f "${output}.partial"
            fi
        else
            echo -e "${YELLOW}  ⚠ Download failed (exit code: $download_status) after ${duration}s${NC}"
            # Don't delete partial file - we'll try to resume
        fi
        
        # Wait before retry (exponential backoff)
        if [ $attempt -lt $max_attempts ]; then
            wait_time=$((2 ** attempt))
            echo -e "${YELLOW}  Waiting ${wait_time}s before retry...${NC}"
            sleep $wait_time
        fi
        
        attempt=$((attempt + 1))
    done
    
    # All attempts failed
    echo -e "${RED}  ✗ Failed to download after $max_attempts attempts${NC}"
    rm -f "${output}.partial"
    return 1
}

# Check if VCFs already exist
EXISTING_COUNT=0
for chr in {1..22} X; do
    VCF_BASE=$(get_vcf_filename "$chr")
    
    # Check for both possible file names (with and without .v2)
    if [ -f "${VCF_BASE}.vcf.gz" ] || [ -f "${VCF_BASE}.v2.vcf.gz" ]; then
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
    if [ "$BUILD" = "grch38" ]; then
        echo -e "${YELLOW}→ Downloading 1000 Genomes High Coverage (GRCh38) VCF files...${NC}"
    else
        echo -e "${YELLOW}→ Downloading 1000 Genomes Phase 3 (GRCh37) VCF files...${NC}"
    fi
    echo "This may take a while (several GB of data)"
    NEED_DOWNLOAD=true
fi

# Download VCFs if needed
if [ "$NEED_DOWNLOAD" = true ]; then
    FAILED_DOWNLOADS=()
    
    for chr in {1..22} X; do
        VCF_BASE=$(get_vcf_filename "$chr")
        VCF_FILE="${VCF_BASE}.vcf.gz"
        VCF_FILE_V2="${VCF_BASE}.v2.vcf.gz"
        TBI_FILE="${VCF_FILE}.tbi"
        
        # Check if file already exists (with or without .v2)
        if [ -f "$VCF_FILE" ] || [ -f "$VCF_FILE_V2" ]; then
            echo -e "${GREEN}✓ chr${chr} already present${NC}"
        else
            echo -e "${BLUE}→ Downloading chromosome ${chr}...${NC}"
            
            # Download VCF file
            if download_file_robust "${FTP_BASE}/${VCF_FILE}" "$VCF_FILE"; then
                # Download index file (with retry)
                echo -e "${BLUE}  → Downloading index...${NC}"
                if ! download_file_robust "${FTP_BASE}/${TBI_FILE}" "$TBI_FILE"; then
                    echo -e "${YELLOW}  ⚠ Warning: Failed to download index for chr${chr}${NC}"
                    echo -e "${YELLOW}  Will try to create index locally later${NC}"
                fi
                echo -e "${GREEN}✓ chr${chr} download complete${NC}"
            else
                echo -e "${RED}✗ Failed to download chr${chr} after all retries${NC}"
                FAILED_DOWNLOADS+=("chr${chr}")
            fi
        fi
    done
    
    # Check if any downloads failed
    if [ ${#FAILED_DOWNLOADS[@]} -gt 0 ]; then
        echo -e "\n${RED}========================================${NC}"
        echo -e "${RED}Download Failed for ${#FAILED_DOWNLOADS[@]} chromosome(s)${NC}"
        echo -e "${RED}========================================${NC}"
        echo -e "${RED}Failed chromosomes: ${FAILED_DOWNLOADS[*]}${NC}"
        echo -e "\n${YELLOW}Suggestions:${NC}"
        echo -e "1. Check your internet connection"
        echo -e "2. Try running the script again (it will resume incomplete downloads)"
        echo -e "3. Manually download the missing files from:"
        echo -e "   ${FTP_BASE}"
        echo -e "4. Check if the FTP server is accessible:"
        echo -e "   ping ftp.1000genomes.ebi.ac.uk"
        exit 1
    fi
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
        VCF_BASE=$(get_vcf_filename "$chr")
        VCF_FILE="${VCF_BASE}.vcf.gz"
        VCF_FILE_V2="${VCF_BASE}.v2.vcf.gz"
        
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
    if [ "$BUILD" = "grch38" ]; then
        echo -e "\n${BLUE}ℹ Step 2a: Concatenating 23 chromosomes (may take 3-7 minutes)${NC}"
    else
        echo -e "\n${BLUE}ℹ Step 2a: Concatenating 23 chromosomes (may take 5-10 minutes)${NC}"
    fi
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
    "$OUTPUT_FILE" \
    "${ORIGINAL_DIR}/${ALLELES_FILE}"

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}✅ Extraction complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Output file: ${OUTPUT_FILE}"
echo -e "\n${BLUE}📁 Cached files in ${DOWNLOAD_DIR}:${NC}"
echo -e "   • all_chromosomes_concatenated.vcf.gz (~35GB)"
echo -e "   • 1000genomes_55aisnps.vcf.gz (~100KB)"
echo -e "\n${BLUE}ℹ Next time you run this script, it will be MUCH faster using the cached files!${NC}"
echo -e "Ready for analysis with: ./run.sh\n"

