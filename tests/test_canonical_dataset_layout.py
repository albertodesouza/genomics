import json

import numpy as np

from genomics.workflows.dataset_builders.non_longevous.build_non_longevous_dataset import _reorganize_output_structure
from genomics.workflows.dataset_builders.non_longevous.dataset_builder import DatasetMetadataBuilder, IndividualDatasetBuilder
from genomics.workflows.dataset_builders.non_longevous.genomic_dataset import GenomicLongevityDataset


def _write_fasta(path, name, sequence):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{name}\n{sequence}\n")


def test_non_longevous_builder_writes_and_reads_canonical_layout(tmp_path):
    dataset_dir = tmp_path / "dataset"
    sample_id = "HG00096"
    target_name = "MC1R"
    individual_base = dataset_dir / "individuals" / sample_id
    legacy_case_dir = individual_base / f"{sample_id}__{target_name}"

    _write_fasta(legacy_case_dir / "ref.window.fa", "chr16:89919736-89919743", "ACGTACGT")
    _write_fasta(legacy_case_dir / f"{sample_id}.H1.window.fixed.fa", "H1", "ACGTTCGT")
    _write_fasta(legacy_case_dir / f"{sample_id}.H2.window.fixed.fa", "H2", "ACGTACGA")
    (legacy_case_dir / "window_metadata.json").write_text(json.dumps({
        "type": "gene",
        "chromosome": "chr16",
        "start": 89919735,
        "end": 89919743,
        "window_size": 8,
        "outputs": ["RNA_SEQ"],
        "ontologies": ["CL:1000458"],
    }))
    pred_dir = legacy_case_dir / "predictions_H1"
    pred_dir.mkdir(parents=True)
    np.savez_compressed(pred_dir / "rna_seq.npz", values=np.array([1.0, 2.0]))

    assert _reorganize_output_structure(individual_base, sample_id) == [target_name]
    assert not (individual_base / "windows" / target_name / "ref.window.fa").exists()
    assert (dataset_dir / "references" / "windows" / target_name / "ref.window.fa").exists()
    assert (dataset_dir / "references" / "windows" / target_name / "window_metadata.json").exists()

    individual_builder = IndividualDatasetBuilder(
        base_dir=dataset_dir,
        sample_id=sample_id,
        sample_info={
            "FamilyID": "FAM1",
            "SampleID": sample_id,
            "Sex": 1,
            "Population": "GBR",
            "Superpopulation": "EUR",
        },
    )
    individual_builder.create_structure()
    individual_builder.add_window(
        target_name=target_name,
        window_type="gene",
        chromosome="chr16",
        start=89919735,
        end=89919743,
        outputs=["RNA_SEQ"],
        ontologies=["CL:1000458"],
    )
    individual_builder.save_metadata()

    metadata_builder = DatasetMetadataBuilder(dataset_dir=dataset_dir, dataset_name="test_dataset", window_size=8)
    metadata_builder.scan_individuals()
    metadata_builder.save_metadata()

    dataset_metadata = json.loads((dataset_dir / "dataset_metadata.json").read_text())
    assert dataset_metadata["layout_version"] == 1
    assert dataset_metadata["genes"] == [target_name]
    assert dataset_metadata["window_catalog"][target_name]["chromosome"] == "chr16"
    assert (dataset_dir / "layout_metadata.json").exists()

    dataset = GenomicLongevityDataset(dataset_dir=dataset_dir, load_predictions=True, load_sequences=True)
    input_data, output_data = dataset[0]
    window_data = input_data["windows"][target_name]

    assert output_data["sample_id"] == sample_id
    assert window_data["ref_sequence"] == "ACGTACGT"
    assert window_data["h1_sequence"] == "ACGTTCGT"
    assert window_data["h2_sequence"] == "ACGTACGA"
    assert window_data["window_metadata"]["chromosome"] == "chr16"
    assert window_data["predictions_h1"]["rna_seq"].tolist() == [1.0, 2.0]
