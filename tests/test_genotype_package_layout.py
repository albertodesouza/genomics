import importlib


def test_genotype_subpackages_import():
    for module_name in [
        "genotype_based_predictor.data.pipeline",
        "genotype_based_predictor.data.processed_dataset",
        "genotype_based_predictor.data.layout",
        "genotype_based_predictor.alignment.dynamic_indel_alignment",
        "genotype_based_predictor.alignment.bcftools_chain_mapper",
        "genotype_based_predictor.apps.dataset_browser",
        "genotype_based_predictor.apps.genomics_workbench",
        "genotype_based_predictor.experiments.train",
        "genotype_based_predictor.experiments.training",
        "genotype_based_predictor.experiments.evaluation",
        "genotype_based_predictor.tools.sync_bcftools_chain_artifacts",
    ]:
        importlib.import_module(module_name)


def test_legacy_module_wrappers_still_export_symbols():
    legacy_dataset = importlib.import_module("genotype_based_predictor.dataset")
    new_dataset = importlib.import_module("genotype_based_predictor.data.processed_dataset")
    assert legacy_dataset.ProcessedGenomicDataset is new_dataset.ProcessedGenomicDataset

    legacy_pipeline = importlib.import_module("genotype_based_predictor.data_pipeline")
    new_pipeline = importlib.import_module("genotype_based_predictor.data.pipeline")
    assert legacy_pipeline.prepare_data is new_pipeline.prepare_data

    legacy_training = importlib.import_module("genotype_based_predictor.training")
    new_training = importlib.import_module("genotype_based_predictor.experiments.training")
    assert legacy_training.Trainer is new_training.Trainer
