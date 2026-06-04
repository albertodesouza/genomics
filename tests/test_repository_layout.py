from pathlib import Path


def test_canonical_config_layout_exists():
    root = Path(__file__).resolve().parents[1]
    expected_dirs = [
        "configs/genomes_analyzer",
        "configs/predictors/genotype_based",
        "configs/predictors/variant_transformer",
        "configs/predictors/snp_ancestry",
        "configs/workflows/non_longevous_dataset",
        "configs/workflows/longevity_dataset",
        "configs/workflows/alphagenome",
        "configs/converters/vcf_to_23andme",
        "configs/legacy/neural_ancestry_predictor_deprecated",
    ]

    for rel_path in expected_dirs:
        assert (root / rel_path).is_dir(), rel_path

    expected_files = [
        "configs/genomes_analyzer/config_human_30x_low_memory.yaml",
        "configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml",
        "configs/predictors/genotype_based/neural_legacy/genes_1000_all.yaml",
        "configs/predictors/variant_transformer/superpopulation.yaml",
        "configs/predictors/snp_ancestry/default.yaml",
        "configs/workflows/non_longevous_dataset/default.yaml",
        "configs/workflows/alphagenome/default.yaml",
        "configs/converters/vcf_to_23andme/default.yaml",
    ]

    for rel_path in expected_files:
        assert (root / rel_path).is_file(), rel_path


def test_script_category_layout_exists():
    root = Path(__file__).resolve().parents[1]
    for rel_path in [
        "scripts/env",
        "scripts/ops",
        "scripts/maintenance",
        "scripts/diagnostics",
        "scripts/experiments",
        "scripts/dev",
    ]:
        assert (root / rel_path).is_dir(), rel_path

    expected_scripts = [
        "scripts/experiments/run_important_experiments_overnight.sh",
        "scripts/env/start_genomics_universal.sh",
        "scripts/maintenance/fix_reference_mismatch.sh",
        "scripts/diagnostics/check_output_structure.py",
        "scripts/dev/demo_csv_feature.sh",
    ]

    for rel_path in expected_scripts:
        assert (root / rel_path).is_file(), rel_path

    removed_wrappers = [
        "scripts/run_important_experiments_overnight.sh",
        "scripts/start_genomics_universal.sh",
        "scripts/fix_reference_mismatch.sh",
        "scripts/check_output_structure.py",
        "scripts/demo_csv_feature.sh",
    ]
    for rel_path in removed_wrappers:
        assert not (root / rel_path).exists(), rel_path


def test_removed_root_legacy_directories_are_absent():
    root = Path(__file__).resolve().parents[1]
    removed_dirs = [
        "build_non_longevous_dataset",
        "neural_module",
        "neural_longevity_dataset",
        "neural_ancestry_predictor_deprecated",
        "genomics_pipeline",
        "genotype_based_predictor",
        "variant_transformer_predictor",
        "snp_ancestry_predictor",
        "vcf_to_23andme",
        "genomes_analyzer_pipeline",
    ]

    for rel_path in removed_dirs:
        assert not (root / rel_path).exists(), rel_path
