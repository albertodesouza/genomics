#!/usr/bin/env python3
"""
build_window_and_predict.py

Pipeline to extract a 1 Mb genomic window around a given gene for a specific
1000 Genomes sample, apply that sample's variants to the reference (per haplotype),
and optionally run AlphaGenome predictions on the personalized sequences.

Requirements (install on your machine where you will run this):
  - samtools
  - bcftools
  - Python packages: pandas, alphagenome (and its deps)
  - hg38 FASTA with index (.fai)
  - 1000G (NYGC 30x) VCF.gz + .tbi for the chromosome containing the gene

Example:
  # List available output types
  python3 build_window_and_predict.py --list-outputs
  
  # List available ontologies (tissues/cells)
  python3 build_window_and_predict.py --list-tissues --api-key YOUR_KEY
  
  # Search for specific ontology (e.g., brain)
  python3 build_window_and_predict.py --list-tissues --filter-tissue brain
  
  # Run predictions for a single ontology
  python3 build_window_and_predict.py \
    --sample HG00096 \
    --gene CYP2B6 \
    --ref-fasta /path/to/hg38.fa \
    --vcf /path/to/1000G_hg38_chr19.vcf.gz \
    --outdir ./out \
    --predict \
    --outputs ATAC \
    --ontology UBERON:0002107
  
  # Run predictions for multiple ontologies
  python3 build_window_and_predict.py \
    --sample HG00096 \
    --gene CYP2B6 \
    --ref-fasta /path/to/hg38.fa \
    --vcf /path/to/1000G_hg38_chr19.vcf.gz \
    --outdir ./out \
    --predict \
    --outputs CAGE \
    --ontology "UBERON:0002107,CL:0002601,UBERON:0000955"

Common tissue CURIEs:
  - UBERON:0002107  (liver)
  - UBERON:0000955  (brain)
  - UBERON:0000948  (heart)
  - UBERON:0002048  (lung)
  - CL:0000182      (hepatocyte)
  - CL:0000540      (neuron)

Notes:
- We auto-detect the 'chr' prefix from the FASTA index (.fai) and coerce the gene
  chromosome accordingly.
- We ensure output sequences match the specified window_size (e.g., 524288 or 1048576 bp)
  by trimming or padding with reference bases as needed (after consensus).
- AlphaGenome predictions are saved as compressed NumPy arrays (.npz) with metadata (.json).
  Each .npz file contains tracks with one value per nucleotide position.

Output Structure:
  outdir/SAMPLE__GENE/
    ├── ref.window.fa                    # Reference sequence for window
    ├── SAMPLE.H1.window.fixed.fa        # Haplotype 1 consensus sequence
    ├── SAMPLE.H2.window.fixed.fa        # Haplotype 2 consensus sequence
    ├── predictions_H1/                  # AlphaGenome predictions for H1
    │   ├── atac.npz                     # ATAC-seq predictions (NumPy arrays)
    │   └── atac_metadata.json           # Track metadata
    └── predictions_H2/                  # AlphaGenome predictions for H2
        ├── atac.npz
        └── atac_metadata.json

Author: ChatGPT (for Alberto)
Last updated: 2025-11-04
"""

import argparse
import os
import sys
import subprocess
import tempfile
import json
from pathlib import Path
from typing import Tuple, Optional, Dict

import pandas as pd
import numpy as np

# AlphaGenome imports
from alphagenome.data import gene_annotation
from alphagenome.models import dna_client
from alphagenome.data import genome

# Note: SEQUENCE_LENGTH_1MB = 1048576 (2^20), not 1000000
# We no longer use a global constant; window_size is passed dynamically


def run(cmd: list, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and stream stderr/stdout upon failure."""
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        print(f"[ERROR] Command failed: {' '.join(cmd)}", file=sys.stderr)
        print(proc.stdout, file=sys.stderr)
        print(proc.stderr, file=sys.stderr)
        raise subprocess.CalledProcessError(proc.returncode, cmd, proc.stdout, proc.stderr)
    return proc


def detect_chr_prefix(fasta_fai: Path) -> str:
    """Detect whether FASTA uses 'chr' prefix (returns 'chr' or '')."""
    with open(fasta_fai, 'r') as f:
        first = f.readline().strip().split('\t')[0]
    return 'chr' if first.startswith('chr') else ''


def coerce_chromosome_name(chrom: str, desired_prefix: str) -> str:
    """Coerce a chromosome name to match FASTA prefix style."""
    has_chr = chrom.startswith('chr')
    if desired_prefix == 'chr' and not has_chr:
        return 'chr' + chrom
    if desired_prefix == '' and has_chr:
        return chrom.replace('chr', '', 1)
    return chrom


def read_chromosome_sizes(fasta_fai: Path) -> Dict[str, int]:
    """Read chromosome sizes from FASTA index (.fai) file."""
    chrom_sizes = {}
    with open(fasta_fai, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_name = parts[0]
                chrom_length = int(parts[1])
                chrom_sizes[chrom_name] = chrom_length
    return chrom_sizes


def load_gtf_feather(gtf_feather: Optional[str]) -> pd.DataFrame:
    """Load GENCODE GTF feather. If none provided, use public URL."""
    default_url = (
        "https://storage.googleapis.com/alphagenome/reference/gencode/"
        "hg38/gencode.v46.annotation.gtf.gz.feather"
    )
    src = gtf_feather or default_url
    print(f"[INFO] Loading GTF feather from: {src}")
    gtf = pd.read_feather(src)
    return gtf


def get_gene_interval(gtf: pd.DataFrame, gene_symbol: Optional[str], gene_id: Optional[str]):
    """Return a genome.Interval for the gene (symbol OR id)."""
    if gene_symbol:
        interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_symbol)
    elif gene_id:
        interval = gene_annotation.get_gene_interval(gtf, gene_id=gene_id)
    else:
        raise ValueError("Provide either --gene or --gene-id")
    return interval


def load_snp_list(snp_file: Path, chr_prefix: str) -> list:
    """
    Load SNP list from file and return list of (chrom, start, end, rsID, alleles).
    Expected format: tab-delimited with header
    ALFRED_UID  dbSNP_rsnumber  chrom  chrom_pos  alleles
    
    Returns:
        List of tuples: (chrom, start, end, rsID, alleles)
    """
    snps = []
    with open(snp_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            alfred_uid, rs_id, chrom, pos, alleles = parts[:5]
            chrom = coerce_chromosome_name(chrom, chr_prefix)
            pos = int(pos)
            # Center 1 Mb window on SNP position
            start = max(0, pos - 500000)
            end = pos + 500000
            snps.append((chrom, start, end, rs_id, alleles))
    return snps


def load_gene_list(gene_file: Path) -> list:
    """
    Load gene list from file (one gene symbol or ID per line).
    
    Supports inline comments:
      - Lines starting with # are ignored (full-line comments)
      - Text after # on a line is ignored (inline comments)
      - Empty lines are ignored
    
    Example file format:
      SLC24A5 # solute carrier family 24 member 5
      SLC45A2 # solute carrier family 45 member 2
      OCA2    # OCA2 melanosomal transmembrane protein
    """
    genes = []
    with open(gene_file, 'r') as f:
        for line in f:
            # Remove inline comments (everything after #)
            line = line.split('#')[0].strip()
            
            # Skip empty lines
            if line:
                genes.append(line)
    return genes


def to_region_1based(chrom_0based: str, start0: int, end0: int, 
                    chrom_size: Optional[int] = None) -> Tuple[str, int, int]:
    """
    Convert 0-based [start, end) to 1-based inclusive region string.
    Clips coordinates to valid range [0, chrom_size).
    
    Args:
        chrom_0based: Chromosome name
        start0: Start position (0-based, may be negative)
        end0: End position (0-based, exclusive)
        chrom_size: Chromosome size for validation (optional)
    
    Returns:
        Tuple of (region_string, actual_start_0based, actual_end_0based)
    """
    # Clip start to >= 0
    actual_start_0based = max(0, start0)
    
    # Clip end to chromosome size if provided
    actual_end_0based = end0
    if chrom_size is not None:
        actual_end_0based = min(end0, chrom_size)
    
    # Convert to 1-based
    region_str = f"{chrom_0based}:{actual_start_0based + 1}-{actual_end_0based}"
    
    return region_str, actual_start_0based, actual_end_0based


def get_gtf_cache_path(outdir: Path, cache_dir: Optional[Path] = None) -> Path:
    """
    Get path to GTF cache file.
    
    Args:
        outdir: Default output directory (used if cache_dir is None)
        cache_dir: Optional dedicated cache directory (for shared cache across samples)
    
    Returns:
        Path to gtf_cache.feather
    """
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir / "gtf_cache.feather"
    return outdir / "gtf_cache.feather"


def load_gtf_with_cache(gtf_feather: Optional[str], cache_path: Path) -> pd.DataFrame:
    """Load GTF with caching to avoid repeated downloads."""
    if cache_path.exists():
        print(f"[INFO] Loading cached GTF from: {cache_path}")
        return pd.read_feather(cache_path)
    
    # No cache, load from source (slow)
    print("[INFO] No GTF cache found. Downloading/loading GTF (this may take a while)...")
    gtf = load_gtf_feather(gtf_feather)
    
    # Save to cache
    print(f"[INFO] Saving GTF cache to: {cache_path}")
    gtf.to_feather(cache_path)
    
    return gtf


def read_fasta_seq_only(path: Path) -> str:
    """Read FASTA (single-record) and return concatenated sequence string (no header)."""
    seq_lines = []
    with open(path, 'r') as f:
        for ln in f:
            if ln.startswith('>'):
                continue
            seq_lines.append(ln.strip())
    return ''.join(seq_lines)


def write_fasta(seq: str, path: Path, header: str = "seq"):
    with open(path, 'w') as f:
        f.write(f">{header}\n")
        # wrap at 60 nt per line
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")


def adjust_to_target_size(consensus_seq: str, ref_window_seq: str, target_size: int,
                         requested_start: int, requested_end: int,
                         actual_start: int, actual_end: int,
                         use_ref_only: bool = False) -> str:
    """
    Ensure final sequence is exactly target_size nt, preserving alignment.
    
    Args:
        consensus_seq: The consensus sequence to adjust (ignored if use_ref_only=True)
        ref_window_seq: Reference window sequence extracted by samtools.
                       When use_ref_only=True, this is the sequence to pad.
                       When use_ref_only=False, this is used as padding source for consensus_seq.
                       Note: ref_window_seq may be shorter than target_size (boundary cases).
        target_size: Target sequence length (e.g., 524288, 1048576)
        requested_start: The start coordinate requested (0-based, may be negative)
        requested_end: The end coordinate requested (0-based)
        actual_start: The actual start coordinate returned by samtools (0-based, >= 0)
        actual_end: The actual end coordinate returned by samtools (0-based)
        use_ref_only: If True, return padded reference; if False, use consensus logic
    
    Returns:
        Adjusted sequence of exactly target_size nucleotides with proper alignment
        
    Note:
        When use_ref_only=True, consensus_seq is completely ignored. The function
        pads ref_window_seq with 'N' at the beginning/end to reach target_size,
        based on the difference between requested and actual coordinates.
    """
    if use_ref_only:
        # Mode: return reference genome only (debug mode)
        ref_len = len(ref_window_seq)
        
        # Calculate N padding needed at beginning (when near chromosome start)
        n_prefix = 0
        if requested_start < 0:
            # requested_start is negative, actual_start will be 0
            n_prefix = actual_start - requested_start  # e.g., 0 - (-193642) = 193642
        
        # Calculate N padding needed at end (when near chromosome end)
        expected_length = requested_end - requested_start
        actual_length = actual_end - actual_start
        n_suffix = 0
        if actual_length < expected_length:
            # Missing bases at the end
            n_suffix = expected_length - actual_length - n_prefix
        
        # Build the properly aligned sequence
        result = 'N' * n_prefix + ref_window_seq + 'N' * n_suffix
        
        # Final validation and adjustment
        if len(result) == target_size:
            return result
        elif len(result) > target_size:
            # Shouldn't happen often, but truncate if needed
            return result[:target_size]
        else:
            # Still missing some bases (edge case)
            still_missing = target_size - len(result)
            return result + 'N' * still_missing
    
    # Original consensus logic (when use_ref_only=False)
    L = len(consensus_seq)
    if L == target_size:
        return consensus_seq
    if L > target_size:
        return consensus_seq[:target_size]
    
    # Need to pad
    needed = target_size - L
    # take from ref window tail following the current length (best-effort alignment)
    pad_from_ref = ref_window_seq[L:L+needed]
    if len(pad_from_ref) < needed:
        pad_from_ref += 'N' * (needed - len(pad_from_ref))
    return consensus_seq + pad_from_ref


def bcftools_consensus_per_hap(ref_window_fa: Path, vcf_window_gz: Path, out_fa: Path, hap: str):
    """Run bcftools consensus for a given hap ('1', '2', or 'IUPAC')."""
    hap_arg = hap if hap in ('1', '2') else 'IUPAC'
    cmd = [
        "bcftools", "consensus",
        "-H", hap_arg,
        "-f", str(ref_window_fa),
        str(vcf_window_gz)
    ]
    proc = run(cmd, check=True)
    # write output (bcftools writes to stdout)
    with open(out_fa, 'w') as f:
        f.write(proc.stdout)


def process_window(
    sample: str,
    target_name: str,
    chrom: str,
    start: int,
    end: int,
    ref_fa: Path,
    vcf_path: Path,
    outdir: Path,
    window_size: int,
    args
) -> Path:
    """
    Process a single genomic window: extract reference, subset VCF, build consensus, run predictions.
    
    Args:
        sample: Sample ID (e.g., HG00096)
        target_name: Name of target (gene symbol, rsID, etc.) for output directory
        chrom: Chromosome name (with correct prefix)
        start: Start position (0-based)
        end: End position (0-based, exclusive)
        ref_fa: Path to reference FASTA
        vcf_path: Path to VCF file
        outdir: Base output directory
        window_size: Expected window size in bp
        args: Parsed arguments from argparse
        
    Returns:
        Path to the case directory
    """
    # 0-based to 1-based region
    # Read chromosome sizes for validation
    ref_fai = Path(str(ref_fa) + ".fai")
    chrom_sizes = read_chromosome_sizes(ref_fai)
    chrom_size = chrom_sizes.get(chrom)
    
    # Clip coordinates and get actual values
    region, actual_start, actual_end = to_region_1based(chrom, start, end, chrom_size)
    
    # Store requested coordinates for later use in adjust_to_target_size
    requested_start = start
    requested_end = end
    
    print(f"\n[INFO] Processing {target_name}: {region} (window={window_size})")
    if requested_start < 0:
        print(f"[INFO] Requested start was negative ({requested_start}), clipped to {actual_start}")
    if chrom_size and requested_end > chrom_size:
        print(f"[INFO] Requested end exceeded chromosome size ({requested_end} > {chrom_size}), clipped to {actual_end}")
    
    # Prepare sample/target-specific subdir
    case_dir = outdir / f"{sample}__{target_name}"
    case_dir.mkdir(parents=True, exist_ok=True)
    
    # 1) Extract reference window
    ref_window_fa = case_dir / "ref.window.fa"
    if ref_window_fa.exists():
        print(f"[INFO] Reference window already exists: {ref_window_fa}")
    else:
        print("[INFO] Extracting reference window with samtools faidx ...")
        proc = run(["samtools", "faidx", str(ref_fa), region], check=True)
        with open(ref_window_fa, "w") as f:
            f.write(proc.stdout)
    
    # 2) Subset VCF to sample+region
    vcf_window = case_dir / f"{sample}.window.vcf.gz"
    vcf_window_tbi = Path(str(vcf_window) + ".tbi")
    if vcf_window.exists() and vcf_window_tbi.exists():
        print(f"[INFO] VCF window already exists: {vcf_window}")
    else:
        print("[INFO] Subsetting VCF to sample+region ...")
        run([
            "bcftools", "view",
            "-s", sample,
            "-r", region,
            "-Oz", "-o", str(vcf_window),
            str(vcf_path)
        ], check=True)
        run(["bcftools", "index", "-t", str(vcf_window)], check=True)
    
    # 2b) Filter unsupported symbolic alleles for consensus
    vcf_cons = case_dir / f"{sample}.window.consensus_ready.vcf.gz"
    vcf_cons_tbi = Path(str(vcf_cons) + ".tbi")
    if vcf_cons.exists() and vcf_cons_tbi.exists():
        print(f"[INFO] Consensus-ready VCF already exists: {vcf_cons}")
    else:
        print("[INFO] Filtering unsupported symbolic ALT alleles for consensus ...")
        run([
            "bcftools", "view",
            "-e", 'ALT~"<" && ALT!="<DEL>" && ALT!="<NON_REF>"',
            "-Oz", "-o", str(vcf_cons),
            str(vcf_window)
        ], check=True)
        run(["bcftools", "index", "-t", str(vcf_cons)], check=True)
    
    # 3) Build consensus sequences
    h1_fa = case_dir / f"{sample}.H1.window.raw.fa"
    h2_fa = case_dir / f"{sample}.H2.window.raw.fa"
    iupac_fa = case_dir / f"{sample}.IUPAC.window.raw.fa"
    
    if h1_fa.exists():
        print(f"[INFO] Consensus H1 already exists: {h1_fa}")
    else:
        print("[INFO] Building consensus H1 ...")
        bcftools_consensus_per_hap(ref_window_fa, vcf_cons, h1_fa, hap='1')
    
    if not args.skip_h2:
        if h2_fa.exists():
            print(f"[INFO] Consensus H2 already exists: {h2_fa}")
        else:
            print("[INFO] Building consensus H2 ...")
            bcftools_consensus_per_hap(ref_window_fa, vcf_cons, h2_fa, hap='2')
    
    if args.also_iupac:
        if iupac_fa.exists():
            print(f"[INFO] Consensus IUPAC already exists: {iupac_fa}")
        else:
            print("[INFO] Building consensus IUPAC ...")
            bcftools_consensus_per_hap(ref_window_fa, vcf_cons, iupac_fa, hap='IUPAC')
    
    # 4) Enforce exact length
    print(f"[INFO] Enforcing exact length ({window_size} bp) for outputs ...")
    ref_seq_raw = read_fasta_seq_only(ref_window_fa)
    
    # Ajustar ref_seq para target_size com padding correto
    # Nota: Passamos ref_seq_raw duas vezes porque quando use_ref_only=True,
    # o primeiro argumento (consensus_seq) é ignorado. A função usa apenas
    # ref_window_seq (segundo argumento) para criar a sequência com padding.
    ref_seq = adjust_to_target_size(
        ref_seq_raw, ref_seq_raw, window_size,
        requested_start, requested_end,
        actual_start, actual_end,
        use_ref_only=True
    )
    
    def fix_and_write(raw_path: Path, fixed_path: Path, header: str, ref_adjusted: str, use_ref: bool = False):
        # Se use_ref=True, não precisamos do arquivo raw
        if use_ref:
            if fixed_path.exists():
                print(f"[INFO] Fixed sequence already exists: {fixed_path}")
                return fixed_path
            # Usar referência já ajustada (não precisa chamar adjust_to_target_size novamente)
            write_fasta(ref_adjusted, fixed_path, header=header)
            print(f"[INFO] Building a debug dataset with reference genome only")
            return fixed_path
        
        # Modo normal: precisa do arquivo raw
        if not raw_path.exists():
            return None
        if fixed_path.exists():
            print(f"[INFO] Fixed sequence already exists: {fixed_path}")
            return fixed_path
        seq = read_fasta_seq_only(raw_path)
        fixed = adjust_to_target_size(
            seq, ref_adjusted, window_size,
            requested_start, requested_end,
            actual_start, actual_end,
            use_ref_only=False
        )
        write_fasta(fixed, fixed_path, header=header)
        return fixed_path
    
    # Aplicar ajuste de tamanho para H1, H2, IUPAC
    # Se args.reference_only, usa apenas referência; senão usa consensus do VCF
    h1_fixed = fix_and_write(h1_fa, case_dir / f"{sample}.H1.window.fixed.fa", 
                             header=f"{sample}_H1_{target_name}", ref_adjusted=ref_seq, use_ref=args.reference_only)
    h2_fixed = None
    if not args.skip_h2:
        h2_fixed = fix_and_write(h2_fa, case_dir / f"{sample}.H2.window.fixed.fa", 
                                 header=f"{sample}_H2_{target_name}", ref_adjusted=ref_seq, use_ref=args.reference_only)
    
    iupac_fixed = None
    if args.also_iupac:
        iupac_fixed = fix_and_write(iupac_fa, case_dir / f"{sample}.IUPAC.window.fixed.fa", 
                                     header=f"{sample}_IUPAC_{target_name}", ref_adjusted=ref_seq, use_ref=args.reference_only)
    
    # 5) Optional: AlphaGenome predictions
    if args.predict:
        print("[INFO] Running AlphaGenome predictions ...")
        api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
        if not api_key:
            raise RuntimeError("AlphaGenome API key not provided. Use --api-key or set ALPHAGENOME_API_KEY env var.")
        client = dna_client.create(api_key)
        
        # Convert string output names to OutputType objects
        output_names = [o.strip() for o in args.outputs.split(",") if o.strip()]
        requested_outputs = []
        
        for out_name in output_names:
            variations = [
                out_name.upper().replace("-", "_").replace(" ", "_"),
                out_name.upper().replace("-", "").replace(" ", ""),
                out_name.replace("-", "_").replace(" ", "_"),
                out_name.replace("-", "").replace(" ", ""),
            ]
            
            found = False
            for attr_name in variations:
                try:
                    output_type = getattr(dna_client.OutputType, attr_name)
                    requested_outputs.append(output_type)
                    print(f"[INFO] Mapped '{out_name}' -> OutputType.{attr_name}")
                    found = True
                    break
                except AttributeError:
                    continue
            
            if not found:
                print(f"[WARN] Output type '{out_name}' not found. Tried: {variations}", file=sys.stderr)
        
        if not requested_outputs:
            print("[ERROR] No valid output types specified.", file=sys.stderr)
            sys.exit(1)
        
        # Parse ontology terms
        ontology_terms = None
        if args.ontology:
            ontology_list = [o.strip() for o in args.ontology.split(',') if o.strip()]
            valid_ontologies = []
            for ont in ontology_list:
                if ':' in ont:
                    valid_ontologies.append(ont)
            if valid_ontologies:
                ontology_terms = valid_ontologies
        
        if not ontology_terms and not args.all_tissues:
            # Show warning about all tissues
            print(f"\n{'='*70}")
            print(f"[WARNING] No --ontology specified! This will use ALL tissues.")
            print(f"{'='*70}\n")
        
        def load_seq(fa_path: Path) -> Optional[str]:
            if fa_path is None:
                return None
            seq = read_fasta_seq_only(fa_path)
            if len(seq) != args.window_size:
                print(f"[WARN] {fa_path.name} length is {len(seq)} (expected {args.window_size}). Skipping.", file=sys.stderr)
                return None
            return seq
        
        h1_seq = load_seq(h1_fixed) if h1_fixed else None
        h2_seq = load_seq(h2_fixed) if h2_fixed else None
        iu_seq = load_seq(iupac_fixed) if iupac_fixed else None
        
        # Run predictions
        def predict_and_touch(tag: str, seq: Optional[str]):
            if seq is None:
                return
            marker = case_dir / f"prediction_{tag}.ok.txt"
            if marker.exists():
                print(f"[INFO] Prediction for {tag} already completed: {marker}")
                return
            
            print(f"[INFO] Running prediction for {tag}...", flush=True)
            
            import time
            start_time = time.time()
            outputs = client.predict_sequence(seq, requested_outputs=requested_outputs, ontology_terms=ontology_terms)
            elapsed = time.time() - start_time
            
            print(f"[INFO] API call completed in {elapsed:.1f} seconds", flush=True)
            
            # Apply rate limiting delay after API call
            if args.api_rate_limit_delay > 0:
                print(f"[INFO] Rate limiting: waiting {args.api_rate_limit_delay}s before next API call...", flush=True)
                time.sleep(args.api_rate_limit_delay)
            
            predictions_dir = case_dir / f"predictions_{tag}"
            predictions_dir.mkdir(exist_ok=True)
            
            saved_outputs = []
            for output_type in requested_outputs:
                attr_name = str(output_type).split('.')[-1].lower()
                
                try:
                    track_data = getattr(outputs, attr_name, None)
                    if track_data is None:
                        continue
                    
                    arrays_dict = {}
                    metadata_dict = {}
                    
                    if hasattr(track_data, 'values'):
                        values = track_data.values
                        arrays_dict['values'] = values
                    
                    if hasattr(track_data, 'metadata'):
                        metadata = track_data.metadata
                        try:
                            if 'DataFrame' in str(type(metadata)):
                                metadata_dict['metadata'] = metadata.to_dict(orient='records')
                            else:
                                metadata_dict['metadata'] = str(metadata)
                        except:
                            metadata_dict['metadata'] = str(metadata)
                    
                    if arrays_dict:
                        npz_file = predictions_dir / f"{attr_name}.npz"
                        np.savez_compressed(npz_file, **arrays_dict)
                        saved_outputs.append(attr_name)
                    
                    if metadata_dict:
                        meta_file = predictions_dir / f"{attr_name}_metadata.json"
                        with open(meta_file, 'w') as f:
                            json.dump(metadata_dict, f, indent=2)
                        
                except Exception as e:
                    print(f"[WARN] Error processing output {attr_name}: {e}", flush=True)
                    continue
            
            marker.write_text(f"AlphaGenome prediction completed.\nOutputs saved: {', '.join(saved_outputs)}\n")
        
        predict_and_touch("H1", h1_seq)
        predict_and_touch("H2", h2_seq)
        predict_and_touch("IUPAC", iu_seq)
    
    print(f"[DONE] Outputs for {target_name} written to: {case_dir}")
    return case_dir


def main():
    ap = argparse.ArgumentParser(description="Build personalized 1Mb window around genes or SNPs for a 1000G sample and (optionally) run AlphaGenome.")
    
    # Mode selection
    ap.add_argument("--mode", choices=['gene', 'snp'], default='gene', 
                    help="Operating mode: 'gene' (center window on gene) or 'snp' (center window on SNP position)")
    
    # Sample and reference
    ap.add_argument("--sample", help="1000 Genomes sample ID (e.g., HG00096)")
    ap.add_argument("--ref-fasta", help="Path to hg38 FASTA")
    
    # Gene mode options
    ap.add_argument("--gene", help="HGNC gene symbol (e.g., CYP2B6) - for gene mode")
    ap.add_argument("--gene-id", help="ENSEMBL gene id (e.g., ENSG...) - for gene mode")
    ap.add_argument("--gene-list-file", help="Path to file with gene list (one gene symbol or ID per line) - for gene mode")
    
    # SNP mode options
    ap.add_argument("--snp-list-file", help="Path to SNP list file (tab-delimited: ALFRED_UID, dbSNP_rsnumber, chrom, chrom_pos, alleles) - for snp mode")
    ap.add_argument("--vcf", help="Path to chromosome-level 1000G VCF.gz (indexed) containing the sample")
    ap.add_argument("--gtf-feather", help="Path to GENCODE v46 GTF feather (optional; will use public URL if omitted)")
    ap.add_argument("--gtf-cache-dir", help="Directory to store shared GTF cache (optional; defaults to --outdir)")
    ap.add_argument("--window-size", type=int, default=1048576, 
                    help="Window size in bp (must be power of 2); default is 1048576 (1MB). Valid: 524288 (512KB), 1048576 (1MB)")
    ap.add_argument("--outdir", default="./out", help="Output directory")
    ap.add_argument("--predict", action="store_true", help="Run AlphaGenome predictions on haplotypes (H1/H2)")
    ap.add_argument("--skip-h2", action="store_true", help="Skip building H2 (haplotype 2)")
    ap.add_argument("--also-iupac", action="store_true", help="Also build an IUPAC-coded window")
    ap.add_argument("--reference-only", action="store_true", 
                    help="Use reference genome only (ignore VCFs, for debug)")
    ap.add_argument("--api-key", help="AlphaGenome API key (or set ALPHAGENOME_API_KEY env var)")
    ap.add_argument("--ontology", "--tissue", dest="ontology", help="Ontology CURIE(s) for tissue/cell type. Single: UBERON:0002107. Multiple (comma-separated): UBERON:0002107,CL:0002601. If not provided, uses all tissues/cells.")
    ap.add_argument("--outputs", default="RNA-seq", help="Comma-separated requested outputs (e.g., RNA-seq,ATAC-seq)")
    ap.add_argument("--list-outputs", action="store_true", help="List available output types and exit")
    ap.add_argument("--list-tissues", action="store_true", help="List available tissue/cell ontologies and exit (requires API key)")
    ap.add_argument("--filter-tissue", help="Filter tissue list by name (case-insensitive, e.g., 'brain' or 'liver'). Use with --list-tissues")
    ap.add_argument("--all-tissues", action="store_true", help="Skip confirmation when requesting all tissues (WARNING: very slow and memory intensive)")
    ap.add_argument("--api-rate-limit-delay", type=float, default=0.0, help="Delay (in seconds) between AlphaGenome API calls to respect usage limits. Supports float values (e.g., 0.5 for 500ms)")
    args = ap.parse_args()
    
    # Handle --list-outputs early
    if args.list_outputs:
        print("Available OutputType attributes in AlphaGenome:")
        for attr in sorted(dir(dna_client.OutputType)):
            if not attr.startswith('_') and attr.isupper():
                print(f"  {attr}")
        sys.exit(0)
    
    # Handle --list-tissues (requires API key)
    if args.list_tissues:
        api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
        if not api_key:
            print("[ERROR] AlphaGenome API key required. Use --api-key or set ALPHAGENOME_API_KEY env var.", file=sys.stderr)
            sys.exit(1)
        
        print("[INFO] Loading tissue metadata from AlphaGenome (this may take a few seconds)...")
        client = dna_client.create(api_key)
        metadata = client.output_metadata(dna_client.Organism.HOMO_SAPIENS).concatenate()
        
        # Get unique tissues/cells
        tissues = metadata[['ontology_curie', 'biosample_name', 'biosample_type']].drop_duplicates()
        
        # Apply filter if specified
        if args.filter_tissue:
            filter_term = args.filter_tissue.lower()
            tissues = tissues[tissues['biosample_name'].str.lower().str.contains(filter_term, na=False)]
            print(f"[INFO] Filtering by: '{args.filter_tissue}'\n")
        
        tissues = tissues.sort_values('biosample_name')
        
        print(f"\n{'='*80}")
        print(f"Available tissues/cells in AlphaGenome ({len(tissues)} total)")
        print(f"{'='*80}\n")
        print(f"{'CURIE':<25} {'Biosample Name':<45} {'Type':<15}")
        print(f"{'-'*25} {'-'*45} {'-'*15}")
        
        for _, row in tissues.iterrows():
            curie = str(row['ontology_curie'])
            name = str(row['biosample_name'])[:44]  # Truncate if too long
            btype = str(row['biosample_type'])
            print(f"{curie:<25} {name:<45} {btype:<15}")
        
        print(f"\n{'='*80}")
        print(f"Usage:")
        print(f"  Single:   --ontology UBERON:0002107")
        print(f"  Multiple: --ontology UBERON:0002107,CL:0002601")
        print(f"  (--tissue also works for backward compatibility)")
        print(f"{'='*80}\n")
        sys.exit(0)
    
    # Validate required arguments (only if not using --list-* options)
    if not args.sample or not args.ref_fasta or not args.vcf:
        print("[ERROR] The following arguments are required: --sample, --ref-fasta, --vcf", file=sys.stderr)
        print("[HINT] Use --list-outputs or --list-tissues to see available options without running the pipeline.", file=sys.stderr)
        sys.exit(1)
    
    # Validate mode-specific arguments
    if args.mode == 'gene':
        if not args.gene and not args.gene_id and not args.gene_list_file:
            print("[ERROR] Gene mode requires one of: --gene, --gene-id, or --gene-list-file", file=sys.stderr)
            sys.exit(1)
    elif args.mode == 'snp':
        if not args.snp_list_file:
            print("[ERROR] SNP mode requires --snp-list-file", file=sys.stderr)
            sys.exit(1)

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ref_fa = Path(args.ref_fasta).resolve()
    ref_fai = Path(str(ref_fa) + ".fai")
    if not ref_fai.exists():
        print(f"[INFO] FASTA index not found; indexing: {ref_fai}")
        run(["samtools", "faidx", str(ref_fa)], check=True)
    
    # Detect chromosome prefix
    prefix = detect_chr_prefix(ref_fai)
    
    sample = args.sample
    vcf_pattern = args.vcf  # May contain {chrom} placeholder
    
    # Collect targets (genes or SNPs) to process
    targets = []  # List of (target_name, chrom, start, end) tuples
    
    if args.mode == 'gene':
        # Load GTF with caching (only needed for gene mode)
        cache_dir = Path(args.gtf_cache_dir) if args.gtf_cache_dir else None
        gtf_cache_path = get_gtf_cache_path(outdir, cache_dir)
        gtf = load_gtf_with_cache(args.gtf_feather, gtf_cache_path)
        
        # Collect genes to process
        genes_to_process = []
        if args.gene_list_file:
            # Load from file
            gene_list_path = Path(args.gene_list_file).resolve()
            genes_to_process = load_gene_list(gene_list_path)
            print(f"[INFO] Loaded {len(genes_to_process)} genes from {gene_list_path}")
        elif args.gene or args.gene_id:
            # Single gene
            genes_to_process = [args.gene or args.gene_id]
        
        # Get intervals for each gene
        for gene in genes_to_process:
            # Determine if it's a symbol or ID
            if gene.startswith('ENSG'):
                interval = get_gene_interval(gtf, gene_symbol=None, gene_id=gene)
                gene_name = gene
            else:
                interval = get_gene_interval(gtf, gene_symbol=gene, gene_id=None)
                gene_name = gene
            
            # Resize to desired window
            interval = interval.resize(args.window_size)
            
            # Coerce chromosome prefix
            chrom = coerce_chromosome_name(interval.chromosome, prefix)
            targets.append((gene_name, chrom, interval.start, interval.end))
        
        print(f"[INFO] Processing {len(targets)} gene(s) in mode: gene")
    
    elif args.mode == 'snp':
        # Load SNPs from file
        snp_file_path = Path(args.snp_list_file).resolve()
        snp_list = load_snp_list(snp_file_path, prefix)
        
        for chrom, start, end, rs_id, alleles in snp_list:
            targets.append((rs_id, chrom, start, end))
        
        print(f"[INFO] Processing {len(targets)} SNP(s) in mode: snp")
    
    # Process each target
    case_dirs = []
    total_targets = len(targets)
    for idx, (target_name, chrom, start, end) in enumerate(targets, 1):
        print(f"\n[INFO] Processing target {idx}/{total_targets}: {target_name} ({chrom}:{start}-{end})")
        
        # Resolve VCF path for this chromosome
        # Replace {chrom} placeholder if present
        vcf_path_str = vcf_pattern.replace('{chrom}', chrom)
        vcf_path = Path(vcf_path_str).resolve()
        
        # Check if VCF exists for this chromosome
        if not vcf_path.exists():
            print(f"[ERROR] VCF not found for chromosome {chrom}: {vcf_path}", file=sys.stderr)
            print(f"[ERROR] Skipping target {target_name}", file=sys.stderr)
            continue
        
        case_dir = process_window(
            sample=sample,
            target_name=target_name,
            chrom=chrom,
            start=start,
            end=end,
            ref_fa=ref_fa,
            vcf_path=vcf_path,
            outdir=outdir,
            window_size=args.window_size,
            args=args
        )
        case_dirs.append(case_dir)
    
    # Summary
    print(f"\n{'='*80}")
    print(f"[SUMMARY] Processed {len(case_dirs)} target(s) for sample {sample}")
    print(f"[SUMMARY] Output directory: {outdir}")
    for i, case_dir in enumerate(case_dirs, 1):
        print(f"[SUMMARY]   {i}. {case_dir.name}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
