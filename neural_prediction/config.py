from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


@dataclass
class DerivedTargetConfig:
    source_field: str
    exclude_unmapped: bool = False
    class_map: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class DatasetInputConfig:
    hf_dataset_path: str
    alphagenome_outputs: List[str]
    haplotype_mode: str
    window_center_size: int
    downsample_factor: int
    genes_to_use: List[str] = field(default_factory=list)
    normalization_method: str = "zscore"
    processed_cache_dir: Optional[str] = None
    loading_strategy: str = "default"


@dataclass
class OutputConfig:
    prediction_target: str
    known_classes: List[str] = field(default_factory=list)
    derived_targets: Dict[str, DerivedTargetConfig] = field(default_factory=dict)


@dataclass
class ModelConfig:
    type: str = "MLP"
    hidden_layers: List[int] = field(default_factory=lambda: [128, 64])
    dropout_rate: float = 0.0
    activation: str = "relu"
    cnn: Dict[str, Any] = field(default_factory=dict)
    cnn2: Dict[str, Any] = field(default_factory=dict)
    sklearn: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TrainingConfig:
    learning_rate: float
    weight_decay: float
    batch_size: int
    num_epochs: int


@dataclass
class CheckpointingConfig:
    save_frequency: int = 0
    save_best: bool = True
    save_best_loss: bool = True
    load_checkpoint: Optional[str] = None


@dataclass
class WandbConfig:
    use_wandb: bool = False
    project_name: Optional[str] = None
    run_name: Optional[str] = None
    log_frequency: int = 100


@dataclass
class DataSplitConfig:
    train_split: float
    val_split: float
    test_split: float
    random_seed: int
    family_split_mode: str = "family_aware"


@dataclass
class PredictionConfig:
    dataset_input: DatasetInputConfig
    output: OutputConfig
    model: ModelConfig
    training: TrainingConfig
    data_split: DataSplitConfig
    checkpointing: CheckpointingConfig = field(default_factory=CheckpointingConfig)
    wandb: WandbConfig = field(default_factory=WandbConfig)
    mode: str = "train"


@dataclass
class ExperimentPaths:
    root: Path
    models_dir: Path
    reports_dir: Path
    config_copy_path: Path


def load_config(config_path: Path) -> PredictionConfig:
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    return parse_config(raw)


def parse_config(raw: Dict[str, Any]) -> PredictionConfig:
    dataset_input = DatasetInputConfig(**raw["dataset_input"])

    derived_targets = {
        key: DerivedTargetConfig(**value)
        for key, value in raw.get("output", {}).get("derived_targets", {}).items()
    }
    output = OutputConfig(
        prediction_target=raw["output"]["prediction_target"],
        known_classes=list(raw.get("output", {}).get("known_classes", [])),
        derived_targets=derived_targets,
    )

    model = ModelConfig(**raw.get("model", {}))
    training = TrainingConfig(**raw["training"])
    data_split = DataSplitConfig(**raw["data_split"])
    checkpointing = CheckpointingConfig(**raw.get("checkpointing", {}))
    wandb = WandbConfig(**raw.get("wandb", {}))

    return PredictionConfig(
        dataset_input=dataset_input,
        output=output,
        model=model,
        training=training,
        data_split=data_split,
        checkpointing=checkpointing,
        wandb=wandb,
        mode=raw.get("mode", "train"),
    )
