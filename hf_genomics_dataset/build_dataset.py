#!/usr/bin/env python4
"""
build_dataset.py

A true build-from-scratch pipeline that constructs a Hugging Face dataset
by pulling directly from the VCF files and AlphaGenome API, bypassing
massive intermediate filesystem layouts.

It reads the same configuration format as build_non_longevous_dataset.py.


Example:

python3 hf_genomics_dataset/build_dataset.py \
  --config build_non_longevous_dataset/configs/gene_1000_test.yaml \
  --output /tmp/test_build   --chunk-size 3

"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Any

from datasets import Dataset, concatenate_datasets, load_from_disk
from tqdm import tqdm

from hf_genomics_dataset.features import gene_window_features, sample_features #Defines the columns
from hf_genomics_dataset.fileio import write_json, stable_json_digest 
from hf_genomics_dataset.models import ConversionOptions #How python will store the data in disk

# Try to import from build_non_longevous_dataset to reuse selection logic
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "build_non_longevous_dataset"))
try:
    from build_non_longevous_dataset import (
        load_config,
        load_metadata_csv,
        select_samples,
        run_build_window_predict,
    )
except ImportError as e:
    print(f"Error importing legacy build tools: {e}")
    sys.exit(1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a direct Hugging Face dataset from scratch."
    )
    parser.add_argument("--config", required=True, help="YAML configuration file.")
    parser.add_argument("--output", required=True, help="Output directory for HF dataset.")
    parser.add_argument("--chunk-size", type=int, default=10, help="Samples per shard.")
    parser.add_argument("--samples", nargs="+", help="Restringir a amostras especificas (ex: HG00096 HG00100)")
    parser.add_argument("--predict", action="store_true", help="Ativar predicoes do AlphaGenome.")
    return parser.parse_args()


# We will need the logic to convert the temporary files into the required dictionary shape.
from hf_genomics_dataset.builder import (
    _build_sample_record,
    _sample_record_to_row,
    _build_gene_window_record,
    _gene_window_record_to_row,
    _flush_rows,
    _merge_shards,
    _build_manifest,
)
from hf_genomics_dataset.models import SourceDataset, BuildSummary

try:
    from build_non_longevous_dataset import (
        parse_fasta_header,
    )
    from dataset_builder import IndividualDatasetBuilder
    from frog_ancestry_parser import FROGAncestryData
except ImportError:
    pass

def load_frog_data(config: Dict, output_dir: Path) -> Optional[FROGAncestryData]:
    # Ensure FROG data is loaded exactly as in legacy
    try:
        likelihood_file = REPO_ROOT / "FROGAncestryCalc/output/whole_1000genomes_55aisnps_likelihood.txt"
        mapping_file = REPO_ROOT / "FROGAncestryCalc/population_mapping_1000genomes.csv"
        if likelihood_file.exists() and mapping_file.exists():
            return FROGAncestryData(likelihood_file, mapping_file)
        else:
            print(f"[WARN] Arquivos do FROGAncestryCalc não encontrados nos caminhos {likelihood_file} e {mapping_file}")
    except Exception as e:
        print(f"Failed to load FROG data: {e}")
    return None

def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output).resolve()
    
    print(f"[INFO] Starting pure build-from-scratch pipeline using {config_path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    config = load_config(config_path)
    csv_path = Path(config["data_sources"]["metadata_csv"])
    if not csv_path.is_absolute():
        csv_path = config_path.parent / csv_path
    df = load_metadata_csv(csv_path)
    #Select the groups of samples that we want
    selected_df = select_samples(df, config)
    

    if args.samples:
        print(f"[INFO] Filtro manual ativado: mantendo apenas {args.samples}")
        selected_df = selected_df[selected_df["SampleID"].isin(args.samples)]
    
    frog_data = load_frog_data(config, output_dir)
    params = config["build_window_params"]
    if args.predict:
        params["predict"] = True
    
    samples_shards_dir = output_dir / ".samples_shards"
    gene_windows_shards_dir = output_dir / ".gene_windows_shards"
    gene_references_shards_dir = output_dir / ".gene_references_shards"
    samples_shards_dir.mkdir(exist_ok=True)
    gene_windows_shards_dir.mkdir(exist_ok=True)
    gene_references_shards_dir.mkdir(exist_ok=True)
    
    sample_buffer = []
    gene_window_buffer = []
    gene_reference_buffer = []
    samples_shard_count = 0
    gene_windows_shard_count = 0
    gene_references_shard_count = 0
    
    gene_names = set()
    total_samples = 0
    total_gene_windows = 0
    total_gene_references = 0
    gene_reference_index = {}
    
    options = ConversionOptions(
        include_sequences=True, 
        chunk_size=args.chunk_size,
        vcf_pattern=config["data_sources"].get("vcf_pattern"),
        include_predictions=params.get("predict", False)
    )
    
    #Main loop
    for idx, row in selected_df.iterrows():
        sample_id = row["SampleID"]
        print(f"\n[INFO] Directly building sample: {sample_id} ({idx + 1}/{len(selected_df)})")
        
        #Create a temporary directory to store the data for the current sample, without filling the hard disk
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            
            #Create the sample info dictionary
            sample_info = {
                "FamilyID": row.get("FamilyID", "0"),
                "SampleID": sample_id,
                "Sex": int(row["Sex"]),
                "Population": row["Population"],
                "Superpopulation": row["Superpopulation"],
            }
            frog_likelihood = None
            frog_pop_names = None
            if frog_data and frog_data.has_sample(sample_id):
                try:
                    frog_info = frog_data.get_individual_data(sample_id, full_frog_vector=True)
                    frog_likelihood = frog_info["likelihood"]
                    frog_pop_names = frog_info["population_names"]
                except Exception:
                    pass
            
            '''
            Create the individual dataset builder, passing the sample info and the frog data
            This creates the data structure of the current sample, creating the needed directories and
            other important structures for a individual.
            '''
            ind_builder = IndividualDatasetBuilder(
                base_dir=temp_dir,
                sample_id=sample_id,
                sample_info=sample_info,
                frog_likelihood=frog_likelihood,
                frog_population_names=frog_pop_names,
            )
            ind_builder.create_structure()
            
            '''
            Run the build_window_and_predict.py script for the current sample
            In summary, the script does this:

            Locates the gene on the DNA map.

            Cuts a large chunk (1Mb) around it.

            Creates a copy of its DNA (with its mutations).

            Predicts with AlphaGenome
            '''
            success, target_name = run_build_window_predict(sample_id, config, temp_dir)
            if not success:
                print(f"[ERROR] Failed building predictions for {sample_id}, skipping")
                continue
                
            windows_dir = temp_dir / "individuals" / sample_id / "windows"
            if windows_dir.exists():
                for window_name in [d.name for d in windows_dir.iterdir() if d.is_dir()]:
                    window_dir = windows_dir / window_name
                    ref_fasta = window_dir / "ref.window.fa"
                    chromosome, start, end = "unknown", 0, params.get("window_size", 1000000)
                    if ref_fasta.exists():
                        parsed_chrom, parsed_start, parsed_end = parse_fasta_header(ref_fasta)
                        if parsed_chrom is not None:
                            chromosome, start, end = parsed_chrom, parsed_start, parsed_end
                    
                    outputs_str = params.get("outputs", "")
                    outputs = [o.strip() for o in outputs_str.split(",") if o.strip()]
                    ontology_str = params.get("ontology", "")
                    ontologies = [o.strip() for o in ontology_str.split(",") if o.strip()]
                    
                    '''
                    Add the window that we created to the individual dataset builder
                    '''
                    ind_builder.add_window(
                        target_name=window_name,
                        window_type=params.get("mode", "gene"),
                        chromosome=chromosome,
                        start=start,
                        end=end,
                        outputs=outputs,
                        ontologies=ontologies,
                    )
            
            ind_builder.save_metadata()
            
            # --- After we did the temporary Legacy way, we transform it in a HF Dict ---

            # Now extract to HF dicts
            source = SourceDataset(name="live_build", path=temp_dir)
            sample_metadata = ind_builder.get_metadata()
            
            
            # Build the sample record (the structure of a sample)
            sample_record = _build_sample_record(source, sample_metadata)
            sample_row = _sample_record_to_row(sample_record)
            
            # Replace real path with placeholder to indicate it wasn't statically sourced
            sample_row["source_path"] = "<live_build>"
            sample_row["source_sample_path"] = f"<live_build>/individuals/{sample_id}"
            
            sample_buffer.append(sample_row)
            total_samples += 1
            
            if args.chunk_size > 0 and len(sample_buffer) >= args.chunk_size:
                samples_shard_count = _flush_rows(sample_buffer, sample_features(), samples_shards_dir, samples_shard_count)
            
            for gene in sample_record.available_genes:
                # Also handle gene reference (deduplicated)
                if gene not in gene_reference_index:
                    from hf_genomics_dataset.builder import _build_gene_reference_record, _gene_reference_record_to_row
                    ref_record = _build_gene_reference_record(source, sample_id, sample_metadata, gene, options)
                    if ref_record:
                        ref_row = _gene_reference_record_to_row(ref_record)
                        gene_reference_buffer.append(ref_row)
                        gene_reference_index[gene] = True
                        total_gene_references += 1
                        
                        if args.chunk_size > 0 and len(gene_reference_buffer) >= args.chunk_size:
                            from hf_genomics_dataset.features import gene_reference_features
                            gene_references_shard_count = _flush_rows(
                                gene_reference_buffer, gene_reference_features(), gene_references_shards_dir, gene_references_shard_count
                            )

                gene_record = _build_gene_window_record(source, sample_id, sample_metadata, gene, options)
                gene_row = _gene_window_record_to_row(gene_record)
                
                gene_row["source_path"] = "<live_build>"
                gene_row["source_sample_path"] = f"<live_build>/individuals/{sample_id}"
                gene_row["source_window_path"] = f"<live_build>/individuals/{sample_id}/windows/{gene}"
                gene_row["source_metadata_path"] = f"<live_build>/individuals/{sample_id}/individual_metadata.json"
                
                gene_window_buffer.append(gene_row)
                gene_names.add(gene)
                total_gene_windows += 1
                
                if args.chunk_size > 0 and len(gene_window_buffer) >= args.chunk_size:
                    gene_windows_shard_count = _flush_rows(
                        gene_window_buffer,
                        gene_window_features(),
                        gene_windows_shards_dir,
                        gene_windows_shard_count,
                    )
    '''
     _flush_rows save temp data in disk using SHARDS

     The script process each sample, one by one, and save the data 
     in disk using SHARDS when the chunk_size is reached
     This is done to avoid memory issues when processing large datasets
    '''
    if sample_buffer:
        _flush_rows(sample_buffer, sample_features(), samples_shards_dir, samples_shard_count)
    if gene_reference_buffer:
        from hf_genomics_dataset.features import gene_reference_features
        _flush_rows(gene_reference_buffer, gene_reference_features(), gene_references_shards_dir, gene_references_shard_count)
    if gene_window_buffer:
        _flush_rows(gene_window_buffer, gene_window_features(), gene_windows_shards_dir, gene_windows_shard_count)

    print("\n[INFO] Merging shards into final dataset...")
    from hf_genomics_dataset.builder import _merge_shards
    _merge_shards(samples_shards_dir, output_dir / "samples")
    _merge_shards(gene_references_shards_dir, output_dir / "gene_reference")
    _merge_shards(gene_windows_shards_dir, output_dir / "gene_windows")
    
    shutil.rmtree(samples_shards_dir)
    shutil.rmtree(gene_references_shards_dir)
    shutil.rmtree(gene_windows_shards_dir)
    
    summary = BuildSummary(
        source_datasets=["live_build"],
        total_samples=total_samples,
        total_gene_windows=total_gene_windows,
        genes=sorted(gene_names),
        duplicate_samples_skipped=0,
        duplicate_gene_windows_skipped=0,
        warnings=[],
    )
    # manifest
    # the _build_manifest expects a list of SourceDataset. We can fake one.
    from hf_genomics_dataset.builder import _build_manifest
    from hf_genomics_dataset.fileio import write_json
    write_json(output_dir / "manifest.json", _build_manifest(summary, [SourceDataset(name="live_build", path=Path("<live_build>"))]))
    
    print(f"\n[DONE] Built pure Hugging Face dataset at {output_dir}")
    print(f"Total samples: {total_samples}")
    print(f"Total gene windows: {total_gene_windows}")
    print(f"Total gene references: {total_gene_references}")

if __name__ == "__main__":
    main()

