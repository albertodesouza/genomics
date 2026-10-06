import importlib


def test_genomics_namespace_exposes_cli_and_workspace():
    cli = importlib.import_module("genomics.cli")
    workspace = importlib.import_module("genomics.workspace")

    assert callable(cli.cli_main)
    assert workspace.DEFAULT_DATASET_ID == "1kg_high_coverage"


def test_genomics_core_compatibility_namespace():
    data_registry = importlib.import_module("genomics.core.data_registry")
    training_utils = importlib.import_module("genomics.core.training_utils")

    assert data_registry.resolve_dataset("1kg_high_coverage").dataset_id == "1kg_high_coverage"
    assert hasattr(training_utils, "EpochTrainer")


def test_variant_transformer_namespace():
    new_config = importlib.import_module("genomics.predictors.variant_transformer.config")
    new_model = importlib.import_module("genomics.predictors.variant_transformer.model")

    assert hasattr(new_config, "ModelConfig")
    assert hasattr(new_model, "VariantTransformerClassifier")


def test_genotype_based_config_and_models_namespace():
    new_config = importlib.import_module("genomics.predictors.genotype_based.config")
    new_models = importlib.import_module("genomics.predictors.genotype_based.models")

    assert hasattr(new_config, "PipelineConfig")
    assert hasattr(new_models, "NNAncestryPredictor")
    assert hasattr(new_models, "CNNAncestryPredictor")
    assert hasattr(new_models, "CNN2AncestryPredictor")


def test_vcf_to_23andme_converter_namespace():
    new_converter = importlib.import_module("genomics.converters.vcf_to_23andme.converter")

    assert hasattr(new_converter, "convert_vcf_to_23andme")
    assert hasattr(new_converter, "vcf_genotype_to_23andme")


def test_snp_ancestry_namespace():
    new_pipeline = importlib.import_module("genomics.predictors.snp_ancestry.pipeline")

    assert hasattr(new_pipeline, "main")
    assert hasattr(new_pipeline, "load_config")


def test_genomes_analyzer_namespace():
    new_pipeline = importlib.import_module("genomics.workflows.genomes_analyzer.pipeline")
    new_cli = importlib.import_module("genomics.workflows.genomes_analyzer.cli")

    assert hasattr(new_pipeline, "main")
    assert hasattr(new_pipeline, "normalize_config_schema")
    assert hasattr(new_cli, "cli_main")


def test_non_longevous_dataset_builder_namespace():
    new_builder = importlib.import_module("genomics.workflows.dataset_builders.non_longevous.build_non_longevous_dataset")
    new_dataset = importlib.import_module("genomics.workflows.dataset_builders.non_longevous.genomic_dataset")

    assert hasattr(new_builder, "main")
    assert hasattr(new_builder, "evaluate_frog_accuracy")
    assert hasattr(new_dataset, "GenomicLongevityDataset")


def test_alphagenome_workflow_namespace():
    new_module = importlib.import_module("genomics.workflows.alphagenome.neural_module")
    new_integration = importlib.import_module("genomics.workflows.alphagenome.neural_integration")

    assert hasattr(new_module, "AlphaGenomeAnalyzer")
    assert hasattr(new_module, "parse_fasta")
    assert hasattr(new_integration, "extract_sequences_from_vcf")


def test_legacy_packages_are_isolated():
    legacy_longevity = importlib.import_module("legacy.neural_longevity_dataset.neural_longevity_dataset")
    legacy_ancestry = importlib.import_module("legacy.neural_ancestry_predictor_deprecated.sklearn_pca_cache")

    assert hasattr(legacy_longevity, "LongevityDatasetBuilder")
    assert hasattr(legacy_ancestry, "run_cli")
