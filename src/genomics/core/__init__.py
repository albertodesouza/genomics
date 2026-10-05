"""Shared infrastructure. Names are imported lazily (PEP 562) so that importing a light module such as
``genomics.core.data_registry`` -- and therefore ``genomics --help`` or ``genomics visualize`` -- does
not pull in torch, scikit-learn or wandb."""
from importlib import import_module
from typing import Any

_EXPORTS = {
    "DatasetRef": "data_registry",
    "resolve_dataset": "data_registry",
    "dataloader_kwargs": "data_loading",
    "make_data_loader_generator": "data_loading",
    "ExperimentRun": "experiment",
    "setup_experiment_run": "experiment",
    "update_manifest": "experiment",
    "classification_metrics": "metrics",
    "save_results_json": "metrics",
    "make_optimizer": "optim",
    "make_optimizer_from_config": "optim",
    "make_torch_generator": "reproducibility",
    "set_random_seeds": "reproducibility",
    "worker_init_fn": "reproducibility",
    "select_device": "run_utils",
    "select_split_loader": "run_utils",
    "training_manifest_fields": "run_utils",
    "SampleRecord": "splitting",
    "SplitSpec": "splitting",
    "build_class_maps": "targets",
    "target_value": "targets",
    "pad_1d": "torch_collate",
    "pad_2d": "torch_collate",
    "move_to_device": "torch_utils",
    "EpochTrainer": "training_utils",
    "append_history_epoch": "training_utils",
    "make_lr_scheduler": "training_utils",
    "new_training_history": "training_utils",
    "step_lr_scheduler": "training_utils",
    "write_training_history": "training_utils",
    "finish_wandb": "wandb_utils",
    "init_wandb_if_enabled": "wandb_utils",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{module}", __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
