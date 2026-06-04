import importlib
from pathlib import Path


def test_genotype_subpackages_import():
    for module_name in [
        "genomics.predictors.genotype_based.utils",
        "genomics.predictors.genotype_based.data.pipeline",
        "genomics.predictors.genotype_based.data.processed_dataset",
        "genomics.predictors.genotype_based.data.layout",
        "genomics.predictors.genotype_based.alignment.dynamic_indel_alignment",
        "genomics.predictors.genotype_based.alignment.bcftools_chain_mapper",
        "genomics.predictors.genotype_based.tools.sync_bcftools_chain_artifacts",
        "genomics.predictors.genotype_based.analysis.plot_sklearn_pca_variance",
        "genomics.predictors.genotype_based.apps.dataset_browser",
        "genomics.predictors.genotype_based.apps.genomics_workbench",
        "genomics.predictors.genotype_based.experiments.train",
        "genomics.predictors.genotype_based.experiments.training",
        "genomics.predictors.genotype_based.experiments.evaluation",
    ]:
        importlib.import_module(module_name)


def test_old_genotype_wrapper_package_removed():
    assert not (Path(__file__).resolve().parents[1] / "genotype_based_predictor" / "config.py").exists()
