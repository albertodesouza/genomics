#!/usr/bin/env python3
"""
Convert SNP coordinates from GRCh37 to GRCh38 using Ensembl REST API

This script reads SNP coordinates from 55_aisnps_alleles.txt (GRCh37/hg19)
and queries the Ensembl REST API to get GRCh38/hg38 coordinates, creating
a new file with the updated coordinates.

Usage:
    python3 convert_grch37_to_grch38.py

Input:  SNPInfo/55_aisnps_alleles.txt (GRCh37)
Output: SNPInfo/55_aisnps_alleles_grch38.txt (GRCh38)
"""

import requests
import time
import sys
import os
from pathlib import Path

# ANSI color codes
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
BLUE = '\033[0;34m'
NC = '\033[0m'  # No Color


def query_ensembl_grch38(rsid):
    """
    Query Ensembl REST API for GRCh38 coordinates of a given rsID.
    
    Args:
        rsid: SNP reference ID (e.g., 'rs1229984')
    
    Returns:
        Tuple of (chromosome, position) or (None, None) if not found
    """
    url = f"https://rest.ensembl.org/variation/human/{rsid}"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        
        if response.status_code == 429:  # Rate limit exceeded
            print(f"{YELLOW}⚠ Rate limit hit, waiting 5 seconds...{NC}")
            time.sleep(5)
            return query_ensembl_grch38(rsid)  # Retry
        
        if response.status_code != 200:
            print(f"{RED}✗ Failed to fetch {rsid}: HTTP {response.status_code}{NC}")
            return None, None
        
        data = response.json()
        
        # Look for GRCh38 mapping in the mappings list
        if 'mappings' in data:
            for mapping in data['mappings']:
                # Check if this is a GRCh38 mapping (assembly_name contains '38')
                if 'assembly_name' in mapping and 'GRCh38' in mapping['assembly_name']:
                    chrom = mapping.get('seq_region_name', '').replace('chr', '')
                    position = mapping.get('start')
                    
                    if chrom and position:
                        return chrom, position
        
        # If no GRCh38 mapping found, check if there's a location
        if 'mappings' not in data or not data['mappings']:
            print(f"{YELLOW}⚠ No mappings found for {rsid}{NC}")
            return None, None
        
        print(f"{YELLOW}⚠ No GRCh38 mapping found for {rsid}{NC}")
        return None, None
        
    except requests.exceptions.Timeout:
        print(f"{RED}✗ Timeout querying {rsid}{NC}")
        return None, None
    except requests.exceptions.RequestException as e:
        print(f"{RED}✗ Error querying {rsid}: {e}{NC}")
        return None, None
    except Exception as e:
        print(f"{RED}✗ Unexpected error for {rsid}: {e}{NC}")
        return None, None


def convert_coordinates(input_file, output_file):
    """
    Convert SNP coordinates from GRCh37 to GRCh38.
    
    Args:
        input_file: Path to input file with GRCh37 coordinates
        output_file: Path to output file for GRCh38 coordinates
    """
    print(f"{BLUE}{'='*60}{NC}")
    print(f"{BLUE}GRCh37 to GRCh38 Coordinate Conversion{NC}")
    print(f"{BLUE}{'='*60}{NC}\n")
    
    if not os.path.exists(input_file):
        print(f"{RED}Error: Input file not found: {input_file}{NC}")
        sys.exit(1)
    
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}\n")
    
    # Read input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 2:
        print(f"{RED}Error: Input file is empty or has no data{NC}")
        sys.exit(1)
    
    # Parse header
    header = lines[0].strip()
    data_lines = lines[1:]
    
    print(f"Processing {len(data_lines)} SNPs...\n")
    
    # Process each SNP
    converted = []
    failed = []
    
    for i, line in enumerate(data_lines, 1):
        parts = line.strip().split('\t')
        
        if len(parts) < 5:
            print(f"{RED}✗ Skipping malformed line {i}{NC}")
            failed.append(parts[1] if len(parts) > 1 else 'unknown')
            continue
        
        alfred_uid = parts[0]
        rsid = parts[1]
        old_chrom = parts[2]
        old_pos = parts[3]
        alleles = parts[4]
        
        print(f"[{i}/{len(data_lines)}] {rsid} (GRCh37: chr{old_chrom}:{old_pos})...", end=' ')
        
        # Query Ensembl for GRCh38 coordinates
        new_chrom, new_pos = query_ensembl_grch38(rsid)
        
        if new_chrom and new_pos:
            print(f"{GREEN}✓ GRCh38: chr{new_chrom}:{new_pos}{NC}")
            converted.append([alfred_uid, rsid, new_chrom, str(new_pos), alleles])
        else:
            print(f"{RED}✗ Failed{NC}")
            failed.append(rsid)
            # Keep original for now with a warning
            converted.append([alfred_uid, rsid, old_chrom, old_pos, alleles])
        
        # Rate limiting: ~15 requests per second = 0.067s per request
        time.sleep(0.07)
    
    # Write output file
    print(f"\n{BLUE}Writing output file...{NC}")
    with open(output_file, 'w') as f:
        f.write(header + '\n')
        for row in converted:
            f.write('\t'.join(row) + '\n')
    
    print(f"{GREEN}✓ Output written to {output_file}{NC}\n")
    
    # Summary
    print(f"{BLUE}{'='*60}{NC}")
    print(f"{BLUE}Summary{NC}")
    print(f"{BLUE}{'='*60}{NC}")
    print(f"Total SNPs: {len(data_lines)}")
    print(f"{GREEN}Successfully converted: {len(data_lines) - len(failed)}{NC}")
    if failed:
        print(f"{RED}Failed conversions: {len(failed)}{NC}")
        print(f"{YELLOW}Failed rsIDs: {', '.join(failed)}{NC}")
        print(f"\n{YELLOW}⚠ Note: Failed SNPs kept original GRCh37 coordinates{NC}")
        print(f"{YELLOW}⚠ You may need to manually verify these positions{NC}")
    else:
        print(f"{GREEN}All SNPs successfully converted!{NC}")
    print()


def main():
    # Determine paths
    script_dir = Path(__file__).parent.parent
    input_file = script_dir / "SNPInfo" / "55_aisnps_alleles.txt"
    output_file = script_dir / "SNPInfo" / "55_aisnps_alleles_grch38.txt"
    
    # Check if output already exists
    if output_file.exists():
        print(f"{YELLOW}⚠ Output file already exists: {output_file}{NC}")
        response = input("Overwrite? (yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("Aborted.")
            sys.exit(0)
        print()
    
    convert_coordinates(str(input_file), str(output_file))


if __name__ == "__main__":
    main()

