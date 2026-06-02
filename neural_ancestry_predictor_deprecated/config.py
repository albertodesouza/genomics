from __future__ import annotations

from pathlib import Path
from typing import Dict

from genomics_pipeline.config_io import load_yaml, write_json
from genomics_pipeline.data_registry import resolve_dataset


def load_config(config_path: Path) -> Dict:
    """Carrega configuração YAML do neural_ancestry_predictor_deprecated."""
    config = load_yaml(config_path)
    dataset_input = config.get("dataset_input")
    if isinstance(dataset_input, dict) and not dataset_input.get("dataset_dir") and dataset_input.get("dataset_id"):
        dataset_input["dataset_dir"] = str(resolve_dataset(str(dataset_input["dataset_id"])).path)
    return config


def save_config(config: Dict, output_path: Path) -> None:
    """Salva configuração como JSON."""
    write_json(output_path, config)


def generate_experiment_name(config: Dict) -> str:
    """Gera nome do experimento baseado nos parâmetros de configuração."""
    model_type = config["model"].get("type", "NN").lower()
    alphagenome_outputs = "_".join(config["dataset_input"]["alphagenome_outputs"])
    haplotype_mode = config["dataset_input"]["haplotype_mode"]
    window_center_size = config["dataset_input"]["window_center_size"]
    normalization_method = config["dataset_input"].get("normalization_method", "zscore")
    hidden_layers = config["model"]["hidden_layers"]
    hidden_layers_str = "L" + "-".join(map(str, hidden_layers))
    activation = config["model"]["activation"]
    dropout_rate = config["model"]["dropout_rate"]
    optimizer = config["training"]["optimizer"]

    if model_type == "cnn":
        cnn_config = config["model"]["cnn"]
        kernel_size = cnn_config["kernel_size"]
        kernel_str = f"k{kernel_size[0]}x{kernel_size[1]}"
        filters_str = f"f{cnn_config['num_filters']}"
        stride = cnn_config["stride"]
        stride_str = f"s{stride[0]}x{stride[1]}" if isinstance(stride, list) else f"s{stride}"
        padding = cnn_config["padding"]
        padding_str = f"p{padding[0]}x{padding[1]}" if isinstance(padding, list) else f"p{padding}"
        pool_size = cnn_config.get("pool_size")
        pool_str = f"pool{pool_size[0]}x{pool_size[1]}_" if pool_size is not None else ""
        return (
            f"cnn_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
            f"{normalization_method}_{kernel_str}_{filters_str}_{stride_str}_{padding_str}_{pool_str}"
            f"{hidden_layers_str}_{activation}_{dropout_rate}_{optimizer}"
        )

    if model_type == "cnn2":
        cnn2_config = config["model"].get("cnn2", {})
        num_filters_s1 = cnn2_config.get("num_filters_stage1", 16)
        kernel_s1 = cnn2_config.get("kernel_stage1", [6, 32])
        stage1_str = f"s1k{kernel_s1[0]}x{kernel_s1[1]}f{num_filters_s1}"
        stage2_str = f"s2f{cnn2_config.get('num_filters_stage2', 32)}"
        stage3_str = f"s3f{cnn2_config.get('num_filters_stage3', 64)}"
        pool_str = f"gp{cnn2_config.get('global_pool_type', 'max')}"
        fc_str = f"fc{cnn2_config.get('fc_hidden_size', 128)}"
        return (
            f"cnn2_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
            f"{normalization_method}_{stage1_str}_{stage2_str}_{stage3_str}_{pool_str}_{fc_str}_"
            f"{hidden_layers_str}_{activation}_{dropout_rate}_{optimizer}"
        )

    if model_type in ("svm", "rf", "xgboost"):
        sk = config["model"].get("sklearn", {})
        pca_k = sk.get("pca_components")
        pca_str = f"pca{pca_k}" if pca_k is not None else "pca_auto"
        if model_type == "svm":
            svm = sk.get("svm", {})
            c_str = str(svm.get("C", 1.0)).replace(".", "p")
            cal = "cal1" if svm.get("calibrate_probabilities", False) else "cal0"
            tag = f"svm_C{c_str}_{cal}"
        elif model_type == "rf":
            rf = sk.get("random_forest", {})
            md = rf.get("max_depth")
            md_str = f"md{md}" if md is not None else "mdNone"
            tag = f"rf_nt{rf.get('n_estimators', 200)}_{md_str}"
        else:
            xgb = sk.get("xgboost", {})
            lr = str(xgb.get("learning_rate", 0.1)).replace(".", "p")
            tag = f"xgb_nt{xgb.get('n_estimators', 200)}_md{xgb.get('max_depth', 6)}_lr{lr}"
        return (
            f"{model_type}_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
            f"{normalization_method}_{pca_str}_{tag}"
        )

    return (
        f"nn_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
        f"{normalization_method}_{hidden_layers_str}_{activation}_{dropout_rate}_{optimizer}"
    )


def generate_dataset_name(config: Dict) -> str:
    """Gera nome único para cache de dataset processado."""
    alphagenome_outputs = "_".join(config["dataset_input"]["alphagenome_outputs"])
    haplotype_mode = config["dataset_input"]["haplotype_mode"]
    window_center_size = config["dataset_input"]["window_center_size"]
    downsample_factor = config["dataset_input"]["downsample_factor"]
    normalization_method = config["dataset_input"].get("normalization_method", "zscore")
    train_split = config["data_split"]["train_split"]
    val_split = config["data_split"]["val_split"]
    test_split = config["data_split"]["test_split"]
    random_seed = config["data_split"]["random_seed"]
    balancing_strategy = config["data_split"].get("balancing_strategy", "stratified")
    balance_str = "strat" if balancing_strategy == "stratified" else "shuf"
    return (
        f"{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
        f"ds{downsample_factor}_{normalization_method}_{balance_str}_"
        f"split{train_split}-{val_split}-{test_split}_seed{random_seed}"
    )


def get_dataset_cache_dir(config: Dict) -> Path:
    base_cache_dir = Path(config["dataset_input"]["processed_cache_dir"])
    return base_cache_dir / "datasets" / generate_dataset_name(config)


def get_results_dir(config: Dict) -> Path:
    return Path(config["dataset_input"].get("results_dir", config["dataset_input"]["processed_cache_dir"]))
