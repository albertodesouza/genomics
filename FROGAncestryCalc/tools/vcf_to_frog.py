#!/usr/bin/env python3
"""
VCF to FROGAncestryCalc Converter
==================================

Converts VCF files (from 1000 Genomes, WGS, or other sources) to the pipe-delimited
format required by FROGAncestryCalc.

Usage:
    python vcf_to_frog.py <input.vcf.gz> <snp_list.txt> <output.txt>

Arguments:
    input.vcf.gz   - VCF file (can be gzipped)
    snp_list.txt   - File with one rsID per line
    output.txt     - Output file in FROGAncestryCalc format

Author: Modified for FROGAncestryCalc pipeline
License: MIT
"""

import sys
import gzip
import os

def open_file(filename):
    """Opens regular or gzipped file automatically"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')

def convert_genotype(gt_field, ref, alt):
    """
    Converts VCF genotype notation (0/1, 0|1) to allele notation (AG, AT, etc.)
    
    Args:
        gt_field: Genotype field from VCF (e.g., '0/1', '1|1')
        ref: Reference allele
        alt: Alternate allele
        
    Returns:
        Two-character genotype string (e.g., 'AG') or 'NN' for missing
    """
    if gt_field in ['./.', '.|.', '.']:
        return 'NN'
    
    # Handle phased (|) or unphased (/) genotypes
    if '|' in gt_field:
        alleles = gt_field.split('|')
    elif '/' in gt_field:
        alleles = gt_field.split('/')
    else:
        return 'NN'
    
    # Convert numeric alleles to actual bases
    result = ''
    for allele in alleles:
        if allele == '0':
            result += ref
        elif allele == '1':
            result += alt.split(',')[0]  # Use first alt allele if multiple
        else:
            return 'NN'
    
    return result

def vcf_to_frog(vcf_file, snp_list_file, output_file):
    """
    Main conversion function
    
    Args:
        vcf_file: Path to VCF file
        snp_list_file: Path to file with SNP IDs (one per line)
        output_file: Path to output file
    """
    
    # Read target SNP list
    print(f"📋 Loading SNP list from {snp_list_file}...")
    with open(snp_list_file) as f:
        target_snps = set(line.strip() for line in f if line.strip())
    
    print(f"   Looking for {len(target_snps)} SNPs")
    
    samples = []
    genotypes = {}  # {sample: {snp: genotype}}
    snp_order = []
    found_snps = set()
    
    print(f"\n🔍 Processing VCF file: {vcf_file}")
    
    with open_file(vcf_file) as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 100000 == 0:
                print(f"   Processed {line_num:,} lines... Found {len(found_snps)} SNPs so far", end='\r')
            
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                # Header with sample names
                fields = line.strip().split('\t')
                samples = fields[9:]  # Samples start at column 10
                print(f"\n✓ Found {len(samples)} samples in VCF")
                
                # Initialize genotype storage
                for sample in samples:
                    genotypes[sample] = {}
                continue
            
            # Variant line
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
                
            chrom, pos, snp_id, ref, alt = fields[0:5]
            
            if snp_id not in target_snps:
                continue
            
            if snp_id in found_snps:
                print(f"\n⚠️  Warning: Duplicate SNP {snp_id} found, keeping first occurrence")
                continue
            
            found_snps.add(snp_id)
            print(f"\n✓ Processing {snp_id} (chr{chrom}:{pos})")
            snp_order.append(snp_id)
            
            # Process genotypes for each sample
            format_field = fields[8].split(':')
            gt_index = format_field.index('GT')
            
            for i, sample in enumerate(samples):
                sample_data = fields[9 + i].split(':')
                gt = sample_data[gt_index]
                genotypes[sample][snp_id] = convert_genotype(gt, ref[0], alt[0])
    
    print(f"\n\n📊 Summary:")
    print(f"   • Found {len(found_snps)} of {len(target_snps)} target SNPs")
    print(f"   • Samples: {len(samples)}")
    
    missing_snps = target_snps - found_snps
    if missing_snps:
        print(f"\n⚠️  Missing SNPs ({len(missing_snps)}):")
        for snp in sorted(missing_snps):
            print(f"      - {snp}")
    
    # Write output file
    print(f"\n💾 Writing output to {output_file}...")
    with open(output_file, 'w') as out:
        # Header
        out.write('Individual|' + '|'.join(snp_order) + '\n')
        
        # Data for each sample
        for sample in samples:
            row = [sample]
            for snp in snp_order:
                row.append(genotypes[sample].get(snp, 'NN'))
            out.write('|'.join(row) + '\n')
    
    print(f"✅ Conversion complete!")
    print(f"   Output file: {output_file}")
    print(f"   Ready for FROGAncestryCalc analysis\n")

def main():
    """Main entry point"""
    if len(sys.argv) != 4:
        print(__doc__)
        print("\nError: Wrong number of arguments")
        print("\nUsage:")
        print("  python vcf_to_frog.py <input.vcf.gz> <snp_list.txt> <output.txt>")
        print("\nExample:")
        print("  python vcf_to_frog.py 1000genomes.vcf.gz tools/aisnps_55_list.txt input/my_samples.txt")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    snp_list_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Validate input files
    if not os.path.exists(vcf_file):
        print(f"❌ Error: VCF file not found: {vcf_file}")
        sys.exit(1)
    
    if not os.path.exists(snp_list_file):
        print(f"❌ Error: SNP list file not found: {snp_list_file}")
        sys.exit(1)
    
    try:
        vcf_to_frog(vcf_file, snp_list_file, output_file)
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()

