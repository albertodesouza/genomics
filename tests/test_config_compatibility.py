from pathlib import Path

import yaml

from genomics.predictors.snp_ancestry.pipeline import load_config as load_snp_config
from genomics.predictors.genotype_based.config import load_config as load_genotype_config


def test_all_canonical_yaml_configs_load():
    for path in sorted(Path("configs").glob("**/*.yaml")):
        with open(path, encoding="utf-8") as f:
            payload = yaml.safe_load(f)
        assert payload is None or isinstance(payload, dict), f"YAML root must be a mapping: {path}"


def test_snp_ancestry_configs_use_canonical_dataset_registry():
    for path in sorted(Path("configs/predictors/snp_ancestry").glob("*.yaml")):
        text = path.read_text(encoding="utf-8")
        assert "/dados/GENOMICS_DATA/top3" not in text

        config = load_snp_config(path)
        input_config = config.get("input", {})
        assert input_config.get("dataset_id") == "1kg_high_coverage"
        assert input_config.get("individuals_dir")
        assert input_config.get("vcf_pattern")

        statistics = config.get("statistics", {})
        prediction = config.get("prediction", {})
        association = config.get("association", {})
        for section in (statistics, prediction, association):
            for key in ("output_dir", "results_dir"):
                value = section.get(key)
                if value:
                    assert not str(value).startswith("/dados/GENOMICS_DATA/top3"), f"{path}:{key} uses legacy top3"


def test_genotype_stability_analysis_config_loads(monkeypatch, tmp_path):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path / "data"))
    config = load_genotype_config(Path("configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml"))

    assert config.stability_analysis.strategy == "repeated_random_split"
    assert config.stability_analysis.n_repeats == 5
    assert config.stability_analysis.stratify is True
