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
- We ensure output sequences are exactly 1,000,000 nt by trimming or padding with
  reference bases as needed (after consensus).
- AlphaGenome predictions are saved as compressed NumPy arrays (.npz) with metadata (.json).
  Each .npz file contains tracks with ~1 million values (one per nucleotide).

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
from typing import Tuple, Optional

import pandas as pd
import numpy as np

# AlphaGenome imports
from alphagenome.data import gene_annotation
from alphagenome.models import dna_client
from alphagenome.data import genome

SEQLEN_1MB = dna_client.SEQUENCE_LENGTH_1MB  # 1_000_000


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


def to_region_1based(chrom_0based: str, start0: int, end0: int) -> str:
    """Convert 0-based [start, end) to 1-based inclusive region string chrom:start-end."""
    return f"{chrom_0based}:{start0 + 1}-{end0}"


def get_gtf_cache_path(outdir: Path) -> Path:
    """Get path to GTF cache file in output directory."""
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


def adjust_to_1mb(consensus_seq: str, ref_window_seq: str) -> str:
    """
    Ensure final sequence is exactly 1,000,000 nt.
    - If longer: trim the end.
    - If shorter: pad with reference bases from the right; if still short, pad 'N'.
    """
    L = len(consensus_seq)
    if L == SEQLEN_1MB:
        return consensus_seq
    if L > SEQLEN_1MB:
        return consensus_seq[:SEQLEN_1MB]

    # Need to pad
    needed = SEQLEN_1MB - L
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


def main():
    ap = argparse.ArgumentParser(description="Build personalized 1Mb window around a gene for a 1000G sample and (optionally) run AlphaGenome.")
    ap.add_argument("--sample", help="1000 Genomes sample ID (e.g., HG00096)")
    ap.add_argument("--gene", help="HGNC gene symbol (e.g., CYP2B6)")
    ap.add_argument("--gene-id", help="ENSEMBL gene id (e.g., ENSG...)")
    ap.add_argument("--ref-fasta", help="Path to hg38 FASTA")
    ap.add_argument("--vcf", help="Path to chromosome-level 1000G VCF.gz (indexed) containing the sample")
    ap.add_argument("--gtf-feather", help="Path to GENCODE v46 GTF feather (optional; will use public URL if omitted)")
    ap.add_argument("--window-size", type=int, default=SEQLEN_1MB, help="Window size; default is AlphaGenome 1Mb")
    ap.add_argument("--outdir", default="./out", help="Output directory")
    ap.add_argument("--predict", action="store_true", help="Run AlphaGenome predictions on haplotypes (H1/H2)")
    ap.add_argument("--skip-h2", action="store_true", help="Skip building H2 (haplotype 2)")
    ap.add_argument("--also-iupac", action="store_true", help="Also build an IUPAC-coded window")
    ap.add_argument("--api-key", help="AlphaGenome API key (or set ALPHAGENOME_API_KEY env var)")
    ap.add_argument("--ontology", "--tissue", dest="ontology", help="Ontology CURIE(s) for tissue/cell type. Single: UBERON:0002107. Multiple (comma-separated): UBERON:0002107,CL:0002601. If not provided, uses all tissues/cells.")
    ap.add_argument("--outputs", default="RNA-seq", help="Comma-separated requested outputs (e.g., RNA-seq,ATAC-seq)")
    ap.add_argument("--list-outputs", action="store_true", help="List available output types and exit")
    ap.add_argument("--list-tissues", action="store_true", help="List available tissue/cell ontologies and exit (requires API key)")
    ap.add_argument("--filter-tissue", help="Filter tissue list by name (case-insensitive, e.g., 'brain' or 'liver'). Use with --list-tissues")
    ap.add_argument("--all-tissues", action="store_true", help="Skip confirmation when requesting all tissues (WARNING: very slow and memory intensive)")
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

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ref_fa = Path(args.ref_fasta).resolve()
    ref_fai = Path(str(ref_fa) + ".fai")
    if not ref_fai.exists():
        print(f"[INFO] FASTA index not found; indexing: {ref_fai}")
        run(["samtools", "faidx", str(ref_fa)], check=True)

    # Load GTF with caching (this is the slow step)
    gtf_cache_path = get_gtf_cache_path(outdir)
    gtf = load_gtf_with_cache(args.gtf_feather, gtf_cache_path)
    
    # Get gene interval and resize to desired window
    interval = get_gene_interval(gtf, gene_symbol=args.gene, gene_id=args.gene_id)
    interval = interval.resize(args.window_size)

    # Coerce chromosome prefix to match FASTA
    prefix = detect_chr_prefix(ref_fai)
    chrom = coerce_chromosome_name(interval.chromosome, prefix)

    # 0-based to 1-based region
    region = to_region_1based(chrom, interval.start, interval.end)
    print(f"[INFO] Region: {region} (window={args.window_size})")

    sample = args.sample
    gene_name = args.gene or args.gene_id or "GENE"

    # Prepare sample/gene-specific subdir
    case_dir = outdir / f"{sample}__{gene_name}"
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
            str(Path(args.vcf).resolve())
        ], check=True)
        run(["bcftools", "index", "-t", str(vcf_window)], check=True)

    # 2b) Filter unsupported symbolic alleles for consensus
    # Exclude any ALT containing "<...>" except "<DEL>" and "<NON_REF>"
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

    # 4) Enforce exact 1 Mb length
    print("[INFO] Enforcing exact length for outputs ...")
    ref_seq = read_fasta_seq_only(ref_window_fa)

    def fix_and_write(raw_path: Path, fixed_path: Path, header: str):
        if not raw_path.exists():
            return None
        if fixed_path.exists():
            print(f"[INFO] Fixed sequence already exists: {fixed_path}")
            return fixed_path
        seq = read_fasta_seq_only(raw_path)
        fixed = adjust_to_1mb(seq, ref_seq)
        write_fasta(fixed, fixed_path, header=header)
        return fixed_path

    h1_fixed = fix_and_write(h1_fa, case_dir / f"{sample}.H1.window.fixed.fa", header=f"{sample}_H1_{gene_name}")
    h2_fixed = None
    if not args.skip_h2:
        h2_fixed = fix_and_write(h2_fa, case_dir / f"{sample}.H2.window.fixed.fa", header=f"{sample}_H2_{gene_name}")

    iupac_fixed = None
    if args.also_iupac:
        iupac_fixed = fix_and_write(iupac_fa, case_dir / f"{sample}.IUPAC.window.fixed.fa", header=f"{sample}_IUPAC_{gene_name}")

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
            # Try multiple naming conventions
            variations = [
                out_name.upper().replace("-", "_").replace(" ", "_"),  # ATAC-seq -> ATAC_SEQ
                out_name.upper().replace("-", "").replace(" ", ""),    # ATAC-seq -> ATACSEQ
                out_name.replace("-", "_").replace(" ", "_"),          # ATAC-seq -> ATAC_seq
                out_name.replace("-", "").replace(" ", ""),            # ATAC-seq -> ATACseq
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
            print("\nAvailable OutputType attributes:", file=sys.stderr)
            for attr in dir(dna_client.OutputType):
                if not attr.startswith('_') and attr.isupper():
                    print(f"  - {attr}", file=sys.stderr)
            sys.exit(1)
        
        # Parse ontology terms (must be CURIEs like "UBERON:0002107" or None for all)
        ontology_terms = None
        if args.ontology:
            # Split by comma to support multiple ontologies
            ontology_list = [o.strip() for o in args.ontology.split(',') if o.strip()]
            
            # Validate each CURIE
            valid_ontologies = []
            for ont in ontology_list:
                if ':' in ont:
                    valid_ontologies.append(ont)
                else:
                    print(f"[WARN] '{ont}' is not a valid CURIE format (expected TYPE:ID). Skipping.")
            
            if valid_ontologies:
                ontology_terms = valid_ontologies
                if len(valid_ontologies) == 1:
                    print(f"[INFO] Using ontology: {valid_ontologies[0]}")
                else:
                    print(f"[INFO] Using {len(valid_ontologies)} ontologies: {', '.join(valid_ontologies)}")
            else:
                print(f"[WARN] No valid CURIE provided.")
                print(f"[INFO] Running predictions for all available tissues/cell types.")
        else:
            print(f"\n{'='*70}")
            print(f"[WARNING] No --ontology specified!")
            print(f"{'='*70}")
            print(f"This will request predictions for ALL available tissues/cell types.")
            print(f"Expected behavior:")
            print(f"  - API call: 5-15 minutes (or more)")
            print(f"  - Memory usage: 10-20 GB")
            print(f"  - Output files: 500 MB - 2 GB per haplotype")
            print(f"  - Hundreds of tracks in the results")
            print(f"\nTo avoid this, specify ontology(ies) with --ontology UBERON:XXXXX")
            print(f"Examples:")
            print(f"  Single:   --ontology UBERON:0002107")
            print(f"  Multiple: --ontology UBERON:0002107,CL:0002601,UBERON:0000955")
            print(f"Common: UBERON:0002107 (liver), UBERON:0000955 (brain), CL:0002601 (senescent cell)")
            print(f"{'='*70}\n")
            
            # Ask for confirmation unless --all-tissues flag is set
            if not args.all_tissues:
                try:
                    response = input("Continue with ALL tissues? [y/N]: ").strip().lower()
                    if response not in ['y', 'yes']:
                        print("[INFO] Aborted by user.")
                        sys.exit(0)
                except (EOFError, KeyboardInterrupt):
                    print("\n[INFO] Aborted by user.")
                    sys.exit(0)
            else:
                print("[INFO] --all-tissues flag set, proceeding without confirmation...")

        def load_seq(fa_path: Path) -> Optional[str]:
            if fa_path is None:
                return None
            seq = read_fasta_seq_only(fa_path)
            if len(seq) != SEQLEN_1MB:
                print(f"[WARN] {fa_path.name} length is {len(seq)} (expected {SEQLEN_1MB}). Skipping.", file=sys.stderr)
                return None
            return seq

        h1_seq = load_seq(h1_fixed) if h1_fixed else None
        h2_seq = load_seq(h2_fixed) if h2_fixed else None
        iu_seq = load_seq(iupac_fixed) if iupac_fixed else None

        # Run predictions and save arrays to disk
        def predict_and_touch(tag: str, seq: Optional[str]):
            if seq is None:
                return
            marker = case_dir / f"prediction_{tag}.ok.txt"
            if marker.exists():
                print(f"[INFO] Prediction for {tag} already completed: {marker}")
                return
            
            print(f"[INFO] Running prediction for {tag}...", flush=True)
            
            # Different messages based on whether ontology_terms is specified
            if ontology_terms is None:
                print(f"[INFO] Calling AlphaGenome API for ALL tissues...", flush=True)
                print(f"[INFO] This may take 5-15 minutes. Please be patient...", flush=True)
            else:
                print(f"[INFO] Calling AlphaGenome API (this may take 30-60 seconds)...", flush=True)
            
            import time
            start_time = time.time()
            outputs = client.predict_sequence(seq, requested_outputs=requested_outputs, ontology_terms=ontology_terms)
            elapsed = time.time() - start_time
            
            print(f"[INFO] API call completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)", flush=True)
            
            # Save prediction arrays for each output type
            predictions_dir = case_dir / f"predictions_{tag}"
            predictions_dir.mkdir(exist_ok=True)
            
            saved_outputs = []
            for output_type in requested_outputs:
                # Get the attribute name (e.g., OutputType.ATAC -> 'atac')
                attr_name = str(output_type).split('.')[-1].lower()
                
                try:
                    # Get the TrackData object
                    track_data = getattr(outputs, attr_name, None)
                    
                    if track_data is None:
                        print(f"[WARN] No data for output type: {attr_name}", flush=True)
                        continue
                    
                    # Extract arrays and metadata from TrackData
                    arrays_dict = {}
                    metadata_dict = {}
                    
                    # Get values (the actual data arrays)
                    if hasattr(track_data, 'values'):
                        values = track_data.values
                        arrays_dict['values'] = values
                        print(f"[INFO] Extracting {attr_name} data: shape={values.shape}", flush=True)
                    
                    # Get number of tracks
                    num_tracks = getattr(track_data, 'num_tracks', 0)
                    metadata_dict['num_tracks'] = num_tracks
                    
                    # Get track names
                    if hasattr(track_data, 'names'):
                        names = track_data.names
                        metadata_dict['track_names'] = names if isinstance(names, list) else list(names)
                    
                    # Get metadata
                    if hasattr(track_data, 'metadata'):
                        metadata = track_data.metadata
                        try:
                            # Try to convert to dict/json-serializable format
                            if hasattr(metadata, 'to_dict'):
                                # For pandas DataFrame, use orient='records' for better JSON structure
                                if 'DataFrame' in str(type(metadata)):
                                    metadata_dict['metadata'] = metadata.to_dict(orient='records')
                                else:
                                    metadata_dict['metadata'] = metadata.to_dict()
                            elif isinstance(metadata, (dict, list)):
                                metadata_dict['metadata'] = metadata
                            else:
                                metadata_dict['metadata'] = str(metadata)
                        except Exception as e:
                            metadata_dict['metadata'] = str(metadata)
                    
                    # Get strand info
                    if hasattr(track_data, 'strands'):
                        strands = track_data.strands
                        metadata_dict['strands'] = strands if isinstance(strands, list) else list(strands)
                    
                    # Get ontology terms
                    if hasattr(track_data, 'ontology_terms'):
                        ont_terms = track_data.ontology_terms
                        # Convert OntologyTerm objects to strings for JSON serialization
                        metadata_dict['ontology_terms'] = [str(term) for term in ont_terms]
                    
                    # Save arrays as compressed npz
                    if arrays_dict:
                        npz_file = predictions_dir / f"{attr_name}.npz"
                        np.savez_compressed(npz_file, **arrays_dict)
                        print(f"[INFO] Saved {len(arrays_dict)} tracks to: {npz_file}", flush=True)
                        saved_outputs.append(attr_name)
                    else:
                        print(f"[WARN] No arrays to save for {attr_name}", flush=True)
                    
                    # Save metadata as JSON
                    if metadata_dict:
                        meta_file = predictions_dir / f"{attr_name}_metadata.json"
                        with open(meta_file, 'w') as f:
                            json.dump(metadata_dict, f, indent=2)
                        print(f"[INFO] Saved metadata to: {meta_file}", flush=True)
                        
                except Exception as e:
                    print(f"[WARN] Error processing output {attr_name}: {e}", flush=True)
                    import traceback
                    traceback.print_exc()
                    continue
            
            # Write completion marker
            marker.write_text(f"AlphaGenome prediction completed.\nOutputs saved: {', '.join(saved_outputs)}\n")
            print(f"[INFO] Prediction for {tag} completed and saved to: {predictions_dir}", flush=True)

        predict_and_touch("H1", h1_seq)
        predict_and_touch("H2", h2_seq)
        predict_and_touch("IUPAC", iu_seq)

        print("[INFO] Predictions completed (markers written).")

    print("\n[DONE] Outputs written to:", case_dir)
    if h1_fixed:
        print(" -", h1_fixed)
    if h2_fixed:
        print(" -", h2_fixed)
    if iupac_fixed:
        print(" -", iupac_fixed)
    print(" -", ref_window_fa)
    print(" -", vcf_window)


if __name__ == "__main__":
    main()
