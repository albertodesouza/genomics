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


def test_root_wrappers_were_removed():
    for module_name in [
        "genotype_based_predictor.dataset",
        "genotype_based_predictor.data_pipeline",
        "genotype_based_predictor.training",
    ]:
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError:
            pass
        else:
            raise AssertionError(f"legacy wrapper should not import: {module_name}")
