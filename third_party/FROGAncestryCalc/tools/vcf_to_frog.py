#!/usr/bin/env python3
"""
VCF to FROGAncestryCalc Converter
==================================

Converts VCF files (from 1000 Genomes, WGS, or other sources) to the pipe-delimited
format required by FROGAncestryCalc.

Usage:
    python vcf_to_frog.py <input.vcf.gz> <snp_list.txt> <output.txt> [alleles_file.txt] [--missing-mode {ref,nn}]

Arguments:
    input.vcf.gz      - VCF file (can be gzipped)
    snp_list.txt      - File with one rsID per line
    output.txt        - Output file in FROGAncestryCalc format
    alleles_file.txt  - Optional: File with SNP coordinates and reference alleles

Options:
    --missing-mode ref  - Missing SNPs become REF/REF when reference alleles are available
    --missing-mode nn   - Missing SNPs become NN even when using alleles_file for coordinate mapping

Author: Modified for FROGAncestryCalc pipeline
License: MIT
"""

import argparse
import sys
import gzip
import os

def open_file(filename):
    """Opens regular or gzipped file automatically"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')

def load_reference_alleles(alleles_file):
    """
    Loads reference alleles from SNPInfo alleles file
    
    Args:
        alleles_file: Path to alleles file (e.g., SNPInfo/55_aisnps_alleles.txt)
        
    Returns:
        Dictionary mapping {rsID: ref_allele} where ref_allele is the first base
        from the alleles column (e.g., 'C' from 'C/T')
    """
    ref_alleles = {}
    
    with open(alleles_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:  # Skip header
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            # Format: ALFRED_UID  dbSNP_rsnumber  chrom  chrom_pos  alleles
            rs_id = parts[1]
            alleles = parts[4]
            
            # Extract first allele (reference) from format like "C/T"
            if '/' in alleles:
                ref_allele = alleles.split('/')[0]
                ref_alleles[rs_id] = ref_allele
    
    return ref_alleles

def load_allowed_alleles(alleles_file):
    """
    Loads allowed alleles from SNPInfo alleles file.
    
    Args:
        alleles_file: Path to alleles file (e.g., SNPInfo/55_aisnps_alleles.txt)
        
    Returns:
        Dictionary mapping {rsID: set(allowed_alleles)}
    """
    allowed_alleles = {}
    
    with open(alleles_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:  # Skip header
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            rs_id = parts[1]
            alleles = parts[4]
            if '/' in alleles:
                allowed_alleles[rs_id] = set(alleles.split('/'))
    
    return allowed_alleles

def complement_base(base):
    """Returns the DNA complement of a single base."""
    complements = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        '.': '.',
    }
    return complements.get(base.upper(), base)

def complement_allele_string(allele_string):
    """Complements a comma-separated allele string."""
    if allele_string in ('', '.'):
        return allele_string
    return ','.join(complement_base(allele) for allele in allele_string.split(','))

def normalize_strand_for_panel(snp_id, ref, alt, allowed_alleles):
    """
    Normalizes REF/ALT alleles to match the strand/orientation expected by the panel.
    
    Args:
        snp_id: SNP rsID
        ref: VCF REF allele
        alt: VCF ALT allele string
        allowed_alleles: Dictionary mapping rsID -> set(allowed alleles)
        
    Returns:
        Tuple: (normalized_ref, normalized_alt, strand_mode)
        strand_mode is one of: 'direct', 'complement', 'unknown'
    """
    expected = allowed_alleles.get(snp_id)
    if not expected:
        return ref, alt, 'unknown'
    
    observed = {ref}
    if alt and alt != '.':
        observed.update(alt.split(','))
    
    if observed.issubset(expected):
        return ref, alt, 'direct'
    
    complemented = {complement_base(base) for base in observed}
    if complemented.issubset(expected):
        return complement_base(ref), complement_allele_string(alt), 'complement'
    
    return ref, alt, 'unknown'

def convert_genotype(gt_field, ref, alt):
    """
    Converts VCF genotype notation (0/1, 0|1) to allele notation (AG, AT, etc.)
    
    Args:
        gt_field: Genotype field from VCF (e.g., '0/1', '1|1')
        ref: Reference allele
        alt: Alternate allele
        
    Returns:
        Two-character genotype string (e.g., 'AG') or 'NN' for missing
        Note: Alleles are sorted alphabetically for consistency
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
    bases = []
    for allele in alleles:
        if allele == '0':
            bases.append(ref)
        elif allele == '1':
            bases.append(alt.split(',')[0])  # Use first alt allele if multiple
        else:
            return 'NN'
    
    # Sort alleles alphabetically for consistency (AG instead of GA, CT instead of TC)
    bases.sort()
    return ''.join(bases)

def vcf_to_frog(vcf_file, snp_list_file, output_file, alleles_file=None, missing_mode='ref'):
    """
    Main conversion function
    
    Args:
        vcf_file: Path to VCF file
        snp_list_file: Path to file with SNP IDs (one per line)
        output_file: Path to output file
        alleles_file: Optional path to alleles file for position mapping and reference alleles
        missing_mode: How to handle missing SNPs. 'ref' uses REF/REF when possible,
            'nn' always writes NN for SNPs absent from the VCF.
    """
    
    # Load reference alleles and position mapping if provided
    ref_alleles_from_file = {}
    allowed_alleles_from_file = {}
    pos_to_rsid = {}  # Maps "chr:pos" to rsID
    
    if alleles_file and os.path.exists(alleles_file):
        print(f"📋 Loading reference alleles and positions from {alleles_file}...")
        ref_alleles_from_file = load_reference_alleles(alleles_file)
        allowed_alleles_from_file = load_allowed_alleles(alleles_file)
        
        # Build position to rsID mapping
        with open(alleles_file, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:  # Skip header
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                # Format: ALFRED_UID  dbSNP_rsnumber  chrom  chrom_pos  alleles
                rs_id = parts[1]
                chrom = parts[2]
                pos = parts[3]
                pos_key = f"{chrom}:{pos}"
                pos_to_rsid[pos_key] = rs_id
        
        print(f"   Loaded {len(ref_alleles_from_file)} reference alleles")
        print(f"   Created position mapping for {len(pos_to_rsid)} SNPs")
        if missing_mode == 'nn':
            print("   Missing SNP mode: NN (alleles file will be used only for coordinate mapping)")
        else:
            print("   Missing SNP mode: REF/REF when reference alleles are available")
    
    # Read target SNP list
    print(f"📋 Loading SNP list from {snp_list_file}...")
    with open(snp_list_file) as f:
        target_snps = set(line.strip() for line in f if line.strip())
    
    print(f"   Looking for {len(target_snps)} SNPs")
    
    samples = []
    genotypes = {}  # {sample: {snp: genotype}}
    snp_order = []
    found_snps = set()
    snp_ref_alleles = {}  # {snp_id: ref_allele} - from VCF
    snp_alt_alleles = {}  # {snp_id: alt_allele} - for logging
    snp_strand_modes = {}  # {snp_id: strand handling mode}
    
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
                
            chrom, pos, vcf_id, ref, alt = fields[0:5]
            
            # Look up rsID by position (chrom:pos)
            # Normalize chromosome name (remove 'chr' prefix if present for matching)
            chrom_normalized = chrom.replace('chr', '') if chrom.startswith('chr') else chrom
            pos_key = f"{chrom_normalized}:{pos}"
            snp_id = pos_to_rsid.get(pos_key)
            
            # If no position mapping, fall back to VCF ID field
            if not snp_id:
                snp_id = vcf_id
            
            # Skip if not in target SNPs
            if snp_id not in target_snps:
                continue
            
            if snp_id in found_snps:
                print(f"\n⚠️  Warning: Duplicate SNP {snp_id} found, keeping first occurrence")
                continue
            
            found_snps.add(snp_id)
            snp_order.append(snp_id)

            normalized_ref, normalized_alt, strand_mode = normalize_strand_for_panel(
                snp_id,
                ref[0],
                alt,
                allowed_alleles_from_file,
            )
            
            # Store reference and alt alleles from VCF
            snp_ref_alleles[snp_id] = normalized_ref
            snp_alt_alleles[snp_id] = normalized_alt[0] if normalized_alt else '.'
            snp_strand_modes[snp_id] = strand_mode
            
            print(f"\n✓ Processing {snp_id} (chr{chrom}:{pos})")
            if strand_mode == 'complement':
                print(f"  ↳ Using strand complement to match panel alleles: REF={ref[0]}/{normalized_ref}")
            elif strand_mode == 'unknown':
                expected = '/'.join(sorted(allowed_alleles_from_file.get(snp_id, [])))
                if expected:
                    print(f"  ↳ Warning: observed alleles do not match panel alleles ({expected}); keeping VCF alleles")
            
            # Process genotypes for each sample
            format_field = fields[8].split(':')
            gt_index = format_field.index('GT')
            
            for i, sample in enumerate(samples):
                sample_data = fields[9 + i].split(':')
                gt = sample_data[gt_index]
                genotypes[sample][snp_id] = convert_genotype(gt, normalized_ref, normalized_alt)
    
    # Process missing SNPs and determine how to handle them
    missing_snps = target_snps - found_snps
    snps_with_ref = set()
    snps_without_ref = set()
    
    print(f"\n\n📊 Summary:")
    print(f"   • Found {len(found_snps)} of {len(target_snps)} target SNPs in VCF")
    print(f"   • Samples: {len(samples)}")
    
    if missing_snps:
        print(f"   • Missing from VCF: {len(missing_snps)} SNPs")
        
        # Determine which missing SNPs have reference alleles
        for snp in missing_snps:
            if missing_mode == 'ref' and snp in ref_alleles_from_file:
                snps_with_ref.add(snp)
            else:
                snps_without_ref.add(snp)
        
        if snps_with_ref:
            print(f"     - Will use REF/REF (homozygous reference): {len(snps_with_ref)} SNPs")
        if snps_without_ref:
            print(f"     - Will use NN (no data): {len(snps_without_ref)} SNPs")
    
    # Detailed processing log for each SNP
    print(f"\n📋 Detailed SNP Processing:")
    print(f"{'='*70}")
    
    # Create complete SNP list (found + missing, in original order from target_snps)
    # We need to maintain order from snp_list_file
    with open(snp_list_file) as f:
        all_snps_ordered = [line.strip() for line in f if line.strip() and line.strip() in target_snps]
    
    for snp in all_snps_ordered:
        if snp in found_snps:
            # Case 1: Found in VCF
            ref = snp_ref_alleles.get(snp, '?')
            alt = snp_alt_alleles.get(snp, '?')
            strand_mode = snp_strand_modes.get(snp, 'direct')
            
            # Count valid genotypes vs NN for this SNP
            valid_count = sum(1 for s in samples if genotypes[s].get(snp, 'NN') != 'NN')
            nn_count = len(samples) - valid_count
            
            status = f"✓ {snp}: Found in VCF (REF={ref}, ALT={alt})"
            if strand_mode == 'complement':
                status += " [strand complemented]"
            elif strand_mode == 'unknown':
                status += " [panel allele mismatch]"
            if nn_count > 0:
                status += f" - {valid_count} valid, {nn_count} NN"
            else:
                status += f" - all {valid_count} samples have data"
            print(status)
            
        elif missing_mode == 'ref' and snp in ref_alleles_from_file:
            # Case 2: Not in VCF, but have reference allele
            ref = ref_alleles_from_file[snp]
            print(f"⚠ {snp}: Not in VCF, using REF/REF ({ref}{ref}) for all {len(samples)} samples")
            snps_with_ref.add(snp)
        else:
            # Case 3: Not in VCF and no reference allele
            if missing_mode == 'nn' and snp in ref_alleles_from_file:
                print(f"⚠ {snp}: Not in VCF, using NN due to --missing-mode nn for all {len(samples)} samples")
            else:
                print(f"✗ {snp}: Not in VCF and no reference available - using NN for all {len(samples)} samples")
            snps_without_ref.add(snp)
    
    print(f"{'='*70}")
    
    # Write output file
    print(f"\n💾 Writing output to {output_file}...")
    with open(output_file, 'w') as out:
        # Header - use complete ordered list
        out.write('Individual|' + '|'.join(all_snps_ordered) + '\n')
        
        # Data for each sample
        for sample in samples:
            row = [sample]
            for snp in all_snps_ordered:
                if snp in found_snps:
                    # Case 1: SNP was in VCF - use stored genotype
                    row.append(genotypes[sample].get(snp, 'NN'))
                elif missing_mode == 'ref' and snp in ref_alleles_from_file:
                    # Case 2: SNP not in VCF but we have reference - use REF/REF
                    ref = ref_alleles_from_file[snp]
                    row.append(ref + ref)
                else:
                    # Case 3: SNP not in VCF and no reference - use NN
                    row.append('NN')
            out.write('|'.join(row) + '\n')
    
    print(f"\n✅ Conversion complete!")
    print(f"   Output file: {output_file}")
    print(f"   Summary:")
    print(f"     - SNPs found in VCF: {len(found_snps)}")
    if snps_with_ref:
        print(f"     - SNPs filled with REF/REF: {len(snps_with_ref)}")
    if snps_without_ref:
        print(f"     - SNPs marked as NN: {len(snps_without_ref)}")
    print(f"   Ready for FROGAncestryCalc analysis\n")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Convert a VCF into FROGAncestryCalc format. When an alleles file is "
            "provided, it is used to map GRCh coordinates to rsIDs; missing SNPs can "
            "be written as REF/REF or NN."
        )
    )
    parser.add_argument("vcf_file", help="VCF file (plain text or gzipped)")
    parser.add_argument("snp_list_file", help="File with one rsID per line")
    parser.add_argument("output_file", help="Output file in FROGAncestryCalc format")
    parser.add_argument(
        "alleles_file",
        nargs="?",
        help="Optional SNPInfo alleles file for coordinate mapping and reference alleles",
    )
    parser.add_argument(
        "--missing-mode",
        choices=["ref", "nn"],
        default="ref",
        help=(
            "How to encode SNPs missing from the VCF. "
            "'ref' uses REF/REF when reference alleles are known; "
            "'nn' always writes NN. Default: ref"
        ),
    )
    return parser.parse_args()


def main():
    """Main entry point"""
    args = parse_args()
    
    vcf_file = args.vcf_file
    snp_list_file = args.snp_list_file
    output_file = args.output_file
    alleles_file = args.alleles_file
    
    # Validate input files
    if not os.path.exists(vcf_file):
        print(f"❌ Error: VCF file not found: {vcf_file}")
        sys.exit(1)
    
    if not os.path.exists(snp_list_file):
        print(f"❌ Error: SNP list file not found: {snp_list_file}")
        sys.exit(1)
    
    if alleles_file and not os.path.exists(alleles_file):
        print(f"❌ Error: Alleles file not found: {alleles_file}")
        sys.exit(1)
    
    try:
        vcf_to_frog(
            vcf_file,
            snp_list_file,
            output_file,
            alleles_file=alleles_file,
            missing_mode=args.missing_mode,
        )
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()

