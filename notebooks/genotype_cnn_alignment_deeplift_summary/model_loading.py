"""Loading CNN2 checkpoints + their test-split data, and running predictions with them."""

from dataclasses import dataclass

import numpy as np
import torch
from sklearn.metrics import classification_report

from genomics.predictors.genotype_based.config import (
    generate_experiment_name,
    get_experiment_runs_dir,
    load_config,
)
from genomics.predictors.genotype_based.data.pipeline import prepare_data
from genomics.predictors.genotype_based.models import CNN2AncestryPredictor

from .logging_utils import quiet_pipeline_logs


def build_model(config, dataset, device):
    input_shape = dataset.get_input_shape()
    num_classes = dataset.get_num_classes()
    if config.model.type.upper() != "CNN2":
        raise ValueError(f"This notebook only compares CNN2 checkpoints, got {config.model.type}")
    return CNN2AncestryPredictor(config, input_shape, num_classes).to(device)


def load_checkpoint(model, checkpoint_path, device):
    checkpoint = torch.load(checkpoint_path, map_location=device)
    state_dict = checkpoint.get("model_state_dict", checkpoint)
    model.load_state_dict(state_dict)
    model.eval()
    return model


def collect_predictions(model, loader, device):
    y_true, y_pred = [], []
    with torch.no_grad():
        for batch in loader:
            features, targets = batch[0], batch[1]
            logits = model(features.to(device))
            preds = logits.argmax(dim=1).cpu().numpy()
            y_true.extend(targets.numpy().tolist())
            y_pred.extend(preds.tolist())
    return np.array(y_true), np.array(y_pred)


@dataclass
class PredictorBundle:
    configs: dict
    experiment_dirs: dict
    datasets: dict
    loaders: dict
    models: dict
    class_names: list


def load_predictor_bundle(config_paths, device, class_names_from):
    """Loads a family of CNN2 checkpoints + their (lazy-loaded, test-split-only) datasets.

    `config_paths` is {name: yaml_path}. `class_names_from` picks which model's dataset supplies
    the class-name labels (any model in the family works -- they all share the same target).
    """
    configs = {name: load_config(path) for name, path in config_paths.items()}
    # Every shipped config uses loading_strategy="preload" (right for training -- avoids
    # per-batch disk I/O), but that means fully materializing train+val+test in RAM. This
    # notebook only ever needs the (much smaller) test split, so force "lazy" loading here --
    # in-memory only, the YAML training configs on disk are untouched.
    for config in configs.values():
        config.data_loading.loading_strategy = "lazy"

    experiment_dirs = {
        name: get_experiment_runs_dir(config) / generate_experiment_name(config)
        for name, config in configs.items()
    }

    for name, config in configs.items():
        di = config.dataset_input
        print(f"{name}: tensor_layout={di.tensor_layout!r} alignment_mapping={di.alignment_mapping!r} "
              f"feature_mode={di.feature_mode!r} indel_mask_positive_value={di.indel_mask_positive_value}")
        print(f"  experiment_dir={experiment_dirs[name]}")
        print(f"  checkpoint exists: {(experiment_dirs[name] / 'models' / 'best_accuracy.pt').exists()}")

    datasets, loaders, models = {}, {}, {}
    # prepare_data's dataset-report step logs cache validation / per-batch progress via a
    # module-level rich Console -- silence it here so 5 models' worth of internal pipeline
    # chatter doesn't bury the one line per model that's actually useful.
    with quiet_pipeline_logs():
        for name, config in configs.items():
            experiment_dir = experiment_dirs[name]
            full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
            datasets[name] = full_ds
            loaders[name] = {"train": train_loader, "val": val_loader, "test": test_loader}

            model = build_model(config, full_ds, device)
            checkpoint_path = experiment_dir / "models" / "best_accuracy.pt"
            models[name] = load_checkpoint(model, checkpoint_path, device)
            n_params = sum(p.numel() for p in model.parameters())
            print(f"{name}: loaded {checkpoint_path.name} ({n_params:,} params), "
                  f"test set size = {len(test_loader.dataset)}")

    class_names = datasets[class_names_from].get_class_names()
    print(f"class_names={class_names}")

    return PredictorBundle(configs, experiment_dirs, datasets, loaders, models, class_names)


def run_predictions(models, loaders, class_names, device):
    """Runs every model in `models` over its own test loader, printing accuracy + a per-class
    classification report, and returns {name: {"y_true":..., "y_pred":...}}."""
    predictions = {}
    for name, model in models.items():
        y_true, y_pred = collect_predictions(model, loaders[name]["test"], device)
        predictions[name] = {"y_true": y_true, "y_pred": y_pred}
        acc = (y_true == y_pred).mean()
        print(f"=== {name} (test, n={len(y_true)}) === accuracy={acc:.4f}")
        print(classification_report(y_true, y_pred, target_names=class_names))
    return predictions
