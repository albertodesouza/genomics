from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, Field, field_validator, model_validator

from genomics_pipeline.config_io import load_yaml, stable_hash
from genomics_pipeline.data_registry import resolve_dataset
from genomics_workspace import DEFAULT_VARIANT_TRANSFORMER_RUNS_ROOT


class DatasetConfig(BaseModel):
    processed_dir: str
    source_dataset_id: Optional[str] = "1kg_high_coverage"
    source_dataset_dir: Optional[str] = None
    results_dir: str = str(DEFAULT_VARIANT_TRANSFORMER_RUNS_ROOT)
    target: str = "superpopulation"
    classes: List[str] = Field(default_factory=lambda: ["AFR", "AMR", "EAS", "EUR", "SAS"])
    max_sequence_length: Optional[int] = None
    truncate_policy: Literal["error", "keep_first", "keep_last"] = "error"
    loading_strategy: Literal["lazy", "preload"] = "lazy"

    @model_validator(mode="before")
    @classmethod
    def resolve_source_dataset_id(cls, values):
        if isinstance(values, dict) and not values.get("source_dataset_dir") and values.get("source_dataset_id"):
            values = dict(values)
            values["source_dataset_dir"] = str(resolve_dataset(str(values["source_dataset_id"])).path)
        return values


class VariantConfig(BaseModel):
    l_max: int = 16
    max_indel_size: int = 50
    unphased_policy: Literal["skip", "error"] = "skip"

    @field_validator("l_max", "max_indel_size")
    @classmethod
    def must_be_positive(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("valor deve ser > 0")
        return value


class ModelConfig(BaseModel):
    d_type: int = 32
    d_hap: int = 16
    d_gene: int = 64
    d_len: int = 16
    d_base: int = 16
    d_allele: int = 128
    d_model: int = 256
    layers: int = 6
    heads: int = 8
    mlp_ratio: int = 4
    dropout: float = 0.1
    rope_base: float = 10000.0

    @model_validator(mode="after")
    def validate_heads(self):
        if self.d_model % self.heads != 0:
            raise ValueError("d_model deve ser divisivel por heads")
        if (self.d_model // self.heads) % 2 != 0:
            raise ValueError("head_dim deve ser par para RoPE")
        return self


class TrainingConfig(BaseModel):
    batch_size: int = 32
    num_epochs: int = 100
    learning_rate: float = 3e-4
    weight_decay: float = 0.01
    optimizer: Literal["adamw", "adam", "sgd"] = "adamw"
    num_workers: int = 0
    early_stopping_patience: Optional[int] = None


class DataSplitConfig(BaseModel):
    train_split: float = 0.7
    val_split: float = 0.15
    test_split: float = 0.15
    random_seed: Optional[int] = 13
    family_split_mode: Literal["family_aware", "ignore"] = "family_aware"


class CheckpointingConfig(BaseModel):
    save_frequency: int = 10
    save_during_training: bool = True
    load_checkpoint: Optional[str] = None


class WandbConfig(BaseModel):
    use_wandb: bool = False
    project_name: str = "variant-transformer"
    run_name: Optional[str] = None


class PipelineConfig(BaseModel):
    dataset: DatasetConfig
    variants: VariantConfig = Field(default_factory=VariantConfig)
    model: ModelConfig = Field(default_factory=ModelConfig)
    training: TrainingConfig = Field(default_factory=TrainingConfig)
    data_split: DataSplitConfig = Field(default_factory=DataSplitConfig)
    checkpointing: CheckpointingConfig = Field(default_factory=CheckpointingConfig)
    wandb: WandbConfig = Field(default_factory=WandbConfig)
    mode: Literal["train", "test"] = "train"
    test_dataset: Literal["train", "val", "test"] = "test"


def load_config(config_path: Path) -> PipelineConfig:
    payload = load_yaml(Path(config_path))
    return PipelineConfig.model_validate(payload)


def save_config(config: PipelineConfig, output_path: Path) -> None:
    with open(output_path, "w") as f:
        f.write(config.model_dump_json(indent=2))


def generate_experiment_name(config: PipelineConfig) -> str:
    payload: Dict[str, object] = {
        "target": config.dataset.target,
        "classes": config.dataset.classes,
        "max_sequence_length": config.dataset.max_sequence_length,
        "model": config.model.model_dump(mode="python"),
        "training": {
            "optimizer": config.training.optimizer,
            "learning_rate": config.training.learning_rate,
            "weight_decay": config.training.weight_decay,
        },
    }
    digest = stable_hash(payload, length=10)
    return f"variant_transformer_{config.dataset.target}_d{config.model.d_model}_l{config.model.layers}_h{config.model.heads}_{digest}"


def get_experiment_dir(config: PipelineConfig) -> Path:
    return Path(config.dataset.results_dir) / generate_experiment_name(config)
