# -*- coding: utf-8 -*-
"""
config.py
=========

Configuração tipada para o pipeline `genotype_based_predictor`.

Usa Pydantic v2 para validação automática e mapeamento YAML → objetos.

Uso rápido
----------
::

    from genotype_based_predictor.config import PipelineConfig, load_config

    config = load_config(Path("configs/pigmentation_binary.yaml"))
    # config é um PipelineConfig com submodelos tipados
    print(config.model.type)            # 'CNN2'
    print(config.training.learning_rate) # 0.001
    print(config.data_split.random_seed) # 13

Funções auxiliares
------------------
load_config              : Carrega YAML e retorna PipelineConfig validado
save_config              : Serializa PipelineConfig como JSON
generate_experiment_name : Gera nome único do experimento
generate_dataset_name    : Gera nome único do dataset (para cache)
get_dataset_cache_dir    : Retorna Path do cache do dataset
"""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


# ===========================================================================
# Submodelos de configuração
# ===========================================================================

class DatasetInputConfig(BaseModel):
    """Parâmetros de entrada dos dados genômicos."""

    dataset_dir: str
    """Caminho absoluto para o diretório do dataset."""

    view_path: Optional[str] = None
    """Arquivo JSON opcional com a definição da view lógica."""

    alphagenome_outputs: List[str]
    """Lista de outputs do AlphaGenome a usar (ex: ['atac', 'rna_seq'])."""

    haplotype_mode: Literal["H1", "H2", "H1+H2"] = "H1+H2"
    """Modo de haplótipo: H1, H2 ou ambos concatenados."""

    tensor_layout: Literal["haplotype_channels"] = "haplotype_channels"
    """Layout canônico do tensor por amostra: (2, 4, L)."""

    window_center_size: int = 32768
    """Número de bases do trecho central de cada janela."""

    downsample_factor: int = 1
    """Fator de downsampling aplicado após extração central (1 = sem downsampling)."""

    genes_to_use: Optional[List[str]] = None
    """Subconjunto de genes a usar. None = todos os genes do dataset."""

    sample_ids: Optional[List[str]] = None
    """Lista explícita de amostras a incluir na análise lógica."""

    sample_ids_path: Optional[str] = None
    """Arquivo texto/JSON com IDs de amostras a incluir."""

    superpopulations_to_use: Optional[List[str]] = None
    """Filtro opcional por superpopulação no dataset lógico."""

    populations_to_use: Optional[List[str]] = None
    """Filtro opcional por população no dataset lógico."""

    gene_order: Optional[List[str]] = None
    """Ordem dos genes no tensor (preenchido automaticamente ao carregar cache)."""

    gene_window_metadata: Optional[Dict[str, Any]] = None
    """Metadados das janelas genômicas (preenchido automaticamente)."""

    normalization_method: Literal["zscore", "minmax_keep_zero", "log"] = "zscore"
    """Método de normalização por track."""

    normalization_value: float = 0.0
    """Valor pré-definido para normalização (legado, ignorado em per-track)."""

    ontology_terms: Optional[List[str]] = None
    """CURIEs de ontologia para filtrar tracks (ex: ['CL:0000236'])."""

    selected_track_index: int = 0
    """Índice da track a usar quando o output possui múltiplas tracks e não há metadados úteis."""

    processed_cache_dir: str
    """Diretório raiz para cache de dados processados e experimentos."""

    keep_split_metadata: bool = False
    """Se True, mantém metadados detalhados por split no cache processado."""

    cache_processed_tensors: bool = True
    """Se False, salva apenas normalização e índices de split, sem gerar train/val/test .pt."""

    normalization_params_path: Optional[str] = None
    """Caminho opcional para normalization_params.json pré-computado e compatível."""

    indel_neutral_value: float = 0.0
    """Valor neutro usado para posições sem sinal real no eixo expandido."""

    indel_include_valid_mask: bool = False
    """Se True, concatena também a máscara de validade além de inserção e deleção."""

    @field_validator("downsample_factor")
    @classmethod
    def downsample_must_be_positive(cls, v: int) -> int:
        if v < 1:
            raise ValueError("downsample_factor deve ser >= 1")
        return v

class DerivedTargetConfig(BaseModel):
    """Configuração de um target derivado a partir de outro campo."""

    source_field: str
    """Campo fonte (ex: 'population')."""

    class_map: Dict[str, List[str]] = Field(default_factory=dict)
    """Mapeamento classe → lista de valores do campo fonte."""

    exclude_unmapped: bool = False
    """Se True, exclui amostras sem mapeamento definido."""


class OutputConfig(BaseModel):
    """Configuração do target de predição."""

    prediction_target: str
    """Target a prever: 'superpopulation', 'population', 'frog_likelihood' ou target derivado."""

    known_classes: Optional[List[str]] = None
    """Lista fixa de classes (evita varredura do dataset para descobrir classes)."""

    derived_targets: Dict[str, DerivedTargetConfig] = Field(default_factory=dict)
    """Targets derivados definidos via mapeamento de campo fonte."""


class CNNConfig(BaseModel):
    """Hiperparâmetros da CNN simples (CNNAncestryPredictor)."""

    kernel_size: List[int] = Field(default=[6, 32])
    num_filters: int = 16
    stride: Union[int, List[int]] = 1
    padding: Union[int, List[int]] = 0
    pool_size: Optional[List[int]] = None


class CNN2Config(BaseModel):
    """Hiperparâmetros da CNN multi-estágio (CNN2AncestryPredictor)."""

    num_filters_stage1: int = 16
    kernel_stage1: List[int] = Field(default=[6, 32])
    num_filters_stage2: int = 32
    kernel_stage2: int = 8
    num_filters_stage3: int = 64
    kernel_stage3: int = 3
    global_pool_type: Literal["max", "avg"] = "max"
    fc_hidden_size: int = 128


class SVMConfig(BaseModel):
    """Hiperparâmetros do SVM baseline."""
    C: float = 1.0
    max_iter: int = 20000
    class_weight: Optional[str] = None
    calibrate_probabilities: bool = False
    calibration_cv: int = 3


class RandomForestConfig(BaseModel):
    """Hiperparâmetros do Random Forest baseline."""
    n_estimators: int = 200
    max_depth: Optional[int] = None
    class_weight: Optional[str] = None
    random_state: int = 42
    n_jobs: int = -1


class XGBoostConfig(BaseModel):
    """Hiperparâmetros do XGBoost baseline."""
    n_estimators: int = 200
    max_depth: int = 6
    learning_rate: float = 0.1
    subsample: float = 0.8
    colsample_bytree: float = 0.8
    tree_method: str = "hist"
    random_state: int = 42
    n_jobs: int = -1
    eval_metric: str = "logloss"


class SklearnConfig(BaseModel):
    """Configuração dos baselines sklearn (PCA + classificador)."""
    pca_components: Optional[int] = None
    use_pca_cache: bool = True
    pca_align_n_train: bool = False
    svm: SVMConfig = Field(default_factory=SVMConfig)
    random_forest: RandomForestConfig = Field(default_factory=RandomForestConfig)
    xgboost: XGBoostConfig = Field(default_factory=XGBoostConfig)


class ModelConfig(BaseModel):
    """Configuração completa do modelo."""

    type: Literal["NN", "CNN", "CNN2", "SVM", "RF", "XGBOOST"] = "NN"
    """Tipo de modelo."""

    hidden_layers: List[int] = Field(default=[128, 64])
    """Tamanhos das camadas ocultas do MLP (NN) ou cabeça FC (CNN/CNN2)."""

    activation: Literal["relu", "tanh", "sigmoid"] = "relu"
    """Função de ativação."""

    dropout_rate: float = 0.0
    """Taxa de dropout aplicada após cada camada oculta."""

    cnn: CNNConfig = Field(default_factory=CNNConfig)
    """Hiperparâmetros da CNN simples."""

    cnn2: CNN2Config = Field(default_factory=CNN2Config)
    """Hiperparâmetros da CNN2 multi-estágio."""

    sklearn: SklearnConfig = Field(default_factory=SklearnConfig)
    """Hiperparâmetros dos baselines sklearn."""

    @field_validator("dropout_rate")
    @classmethod
    def dropout_in_range(cls, v: float) -> float:
        if not 0.0 <= v < 1.0:
            raise ValueError("dropout_rate deve estar em [0, 1)")
        return v


class LRSchedulerConfig(BaseModel):
    """Configuração do learning rate scheduler."""

    enabled: bool = False
    type: Literal[
        "plateau", "step", "cosine", "exponential",
        "multistep", "cosine_warm_restarts"
    ] = "plateau"
    mode: Literal["min", "max"] = "min"
    factor: float = 0.5
    patience: int = 10
    min_lr: float = 1e-6
    step_size: int = 30
    gamma: float = 0.1
    milestones: List[int] = Field(default=[30, 60, 90])
    T_max: int = 100
    eta_min: float = 1e-6
    T_0: int = 50
    T_mult: int = 1


class TrainingConfig(BaseModel):
    """Hiperparâmetros de treinamento."""

    optimizer: Literal["adam", "adamw", "sgd"] = "adam"
    learning_rate: float = 0.001
    weight_decay: float = 0.0
    loss_function: Literal["cross_entropy", "mse"] = "cross_entropy"
    batch_size: int = 128
    num_epochs: int = 100
    validation_frequency: int = 1
    early_stopping_patience: Optional[int] = None
    lr_scheduler: LRSchedulerConfig = Field(default_factory=LRSchedulerConfig)

    @field_validator("learning_rate", "weight_decay")
    @classmethod
    def must_be_non_negative(cls, v: float) -> float:
        if v < 0:
            raise ValueError("Deve ser >= 0")
        return v


class DataSplitConfig(BaseModel):
    """Configuração da divisão treino/validação/teste."""

    train_split: float = 0.7
    val_split: float = 0.15
    test_split: float = 0.15
    random_seed: Optional[int] = 42
    family_split_mode: Literal["family_aware", "ignore"] = "family_aware"
    balancing_strategy: Literal["stratified", "shuffle"] = "stratified"
    strict_determinism: bool = True

    @model_validator(mode="after")
    def splits_must_sum_to_one(self) -> "DataSplitConfig":
        total = self.train_split + self.val_split + self.test_split
        if abs(total - 1.0) > 1e-6:
            raise ValueError(
                f"train_split + val_split + test_split deve ser 1.0 (atual: {total:.4f})"
            )
        return self


class CheckpointingConfig(BaseModel):
    """Configuração de checkpoints do modelo."""

    checkpoint_dir: str = "models"
    save_frequency: int = 10
    save_during_training: bool = True
    load_checkpoint: Optional[str] = None


class WandbConfig(BaseModel):
    """Configuração do Weights & Biases."""

    use_wandb: bool = False
    project_name: str = "genotype-based-predictor"
    run_name: Optional[str] = None
    log_frequency: int = 1


class DebugConfig(BaseModel):
    """Configuração de debug e visualização."""

    enable_visualization: bool = False
    force_pca_cache_rebuild: bool = False
    max_samples_per_epoch: Optional[int] = None
    taint_at_cache_save: bool = False
    taint_at_runtime: bool = False
    taint_type: str = "additive"
    taint_value: float = 1.0
    taint_horizontal_size: int = 6
    taint_vertical_size: int = 6
    taint_horizontal_step: int = 100
    taint_vertical_step: int = 6
    visualization: Dict[str, Any] = Field(default_factory=dict)


class DataLoadingConfig(BaseModel):
    """Configuração de estratégia de carregamento dos dados."""

    loading_strategy: Literal["preload", "lazy"] = "preload"
    cache_size: int = 100


# ===========================================================================
# Modelo raiz
# ===========================================================================

class PipelineConfig(BaseModel):
    """
    Configuração completa do pipeline de predição genômica.

    Mapeado automaticamente a partir do YAML via ``load_config()``.

    Exemplo
    -------
    ::

        config = load_config(Path("configs/pigmentation_binary.yaml"))
        print(config.model.type)              # 'CNN2'
        print(config.training.num_epochs)     # 400
        print(config.data_split.random_seed)  # 13
    """

    dataset_input: DatasetInputConfig
    output: OutputConfig
    model: ModelConfig = Field(default_factory=ModelConfig)
    training: TrainingConfig = Field(default_factory=TrainingConfig)
    data_split: DataSplitConfig = Field(default_factory=DataSplitConfig)
    checkpointing: CheckpointingConfig = Field(default_factory=CheckpointingConfig)
    wandb: WandbConfig = Field(default_factory=WandbConfig)
    debug: DebugConfig = Field(default_factory=DebugConfig)
    data_loading: DataLoadingConfig = Field(default_factory=DataLoadingConfig)

    mode: Literal["train", "test"] = "train"
    test_dataset: Literal["train", "val", "test"] = "test"

    class Config:
        # Permite campos extras no YAML sem falhar (ex: chaves legadas)
        extra = "ignore"


# ===========================================================================
# Funções de I/O
# ===========================================================================

def load_config(config_path: Path) -> PipelineConfig:
    """
    Carrega um arquivo YAML e retorna um ``PipelineConfig`` validado.

    A validação Pydantic garante que todos os campos obrigatórios estão
    presentes e que os valores são do tipo correto. Erros de configuração
    geram ``pydantic.ValidationError`` com mensagens claras.

    Parameters
    ----------
    config_path : Path
        Caminho para o arquivo .yaml.

    Returns
    -------
    PipelineConfig
        Objeto de configuração tipado e validado.

    Raises
    ------
    pydantic.ValidationError
        Se o YAML estiver incompleto ou com tipos inválidos.
    FileNotFoundError
        Se o arquivo não existir.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config não encontrado: {config_path}")

    with open(config_path, "r") as f:
        raw = yaml.safe_load(f)

    dataset_input = raw.get("dataset_input", {}) or {}
    view_path = dataset_input.get("view_path")
    if view_path:
        view_file = Path(view_path)
        if not view_file.is_absolute():
            view_file = (config_path.parent / view_file).resolve()
        if not view_file.exists():
            raise FileNotFoundError(f"View não encontrada: {view_file}")
        with open(view_file) as f:
            view_payload = json.load(f)

        merged_dataset_input = dict(view_payload)
        merged_dataset_input.update(dataset_input)
        merged_dataset_input.pop("name", None)
        merged_dataset_input.pop("description", None)
        merged_dataset_input.pop("prediction_target", None)
        merged_dataset_input.pop("known_classes", None)
        merged_dataset_input.pop("derived_targets", None)
        merged_dataset_input.pop("model", None)
        merged_dataset_input.pop("training", None)
        raw["dataset_input"] = merged_dataset_input

    return PipelineConfig.model_validate(raw)


def save_config(config: PipelineConfig, output_path: Path) -> None:
    """
    Serializa um ``PipelineConfig`` como JSON.

    Parameters
    ----------
    config : PipelineConfig
        Objeto de configuração a serializar.
    output_path : Path
        Caminho de destino para o arquivo JSON.
    """
    output_path = Path(output_path)
    with open(output_path, "w") as f:
        f.write(config.model_dump_json(indent=2))


# ===========================================================================
# Geração de nomes (usados para cache e logging)
# ===========================================================================

def generate_dataset_name(config: PipelineConfig) -> str:
    """
    Gera nome único para o dataset baseado nos parâmetros de processamento.

    Experimentos com os mesmos parâmetros de dataset compartilham o mesmo
    cache, evitando reprocessamento desnecessário.

    Returns
    -------
    str
        Nome único do dataset (ex: ``rna_seq_H1+H2_32768_ds1_log_strat_split0.7-0.15-0.15_seed13``).
    """
    outputs = "_".join(config.dataset_input.alphagenome_outputs)
    hap = config.dataset_input.haplotype_mode
    wcs = config.dataset_input.window_center_size
    ds = config.dataset_input.downsample_factor
    norm = config.dataset_input.normalization_method
    tr = config.data_split.train_split
    va = config.data_split.val_split
    te = config.data_split.test_split
    seed = config.data_split.random_seed
    bal = "strat" if config.data_split.balancing_strategy == "stratified" else "shuf"
    di = config.dataset_input
    view_payload = {
        "dataset_dir": str(Path(di.dataset_dir).resolve()),
        "sample_ids": di.sample_ids,
        "sample_ids_path": di.sample_ids_path,
        "superpopulations_to_use": di.superpopulations_to_use,
        "populations_to_use": di.populations_to_use,
        "genes_to_use": di.genes_to_use,
        "alphagenome_outputs": di.alphagenome_outputs,
        "haplotype_mode": di.haplotype_mode,
        "tensor_layout": di.tensor_layout,
        "window_center_size": di.window_center_size,
        "downsample_factor": di.downsample_factor,
        "normalization_method": di.normalization_method,
        "indel_include_valid_mask": di.indel_include_valid_mask,
        "prediction_target": config.output.prediction_target,
        "train_split": tr,
        "val_split": va,
        "test_split": te,
        "random_seed": seed,
        "balancing_strategy": config.data_split.balancing_strategy,
        "family_split_mode": config.data_split.family_split_mode,
    }
    view_hash = hashlib.sha1(json.dumps(view_payload, sort_keys=True).encode("utf-8")).hexdigest()[:12]
    return f"{outputs}_{hap}_{wcs}_ds{ds}_{norm}_{bal}_view{view_hash}"


def get_dataset_cache_dir(config: PipelineConfig) -> Path:
    """
    Retorna o diretório de cache do dataset.

    Formato: ``<processed_cache_dir>/datasets/<dataset_name>``

    Parameters
    ----------
    config : PipelineConfig

    Returns
    -------
    Path
    """
    base = Path(config.dataset_input.processed_cache_dir)
    return base / "datasets" / generate_dataset_name(config)


def generate_experiment_name(config: PipelineConfig) -> str:
    """
    Gera nome único para o experimento baseado nos parâmetros do modelo.

    O nome codifica tipo de modelo, outputs, janela, normalização e
    principais hiperparâmetros para identificação sem abrir arquivos.

    Returns
    -------
    str
        Nome único do experimento.
    """
    m = config.model
    di = config.dataset_input
    model_type = m.type.lower()

    outputs = "_".join(di.alphagenome_outputs)
    hap = di.haplotype_mode
    tensor_layout = di.tensor_layout
    wcs = di.window_center_size
    norm = di.normalization_method
    hl = "L" + "-".join(str(h) for h in m.hidden_layers)
    act = m.activation
    dr = m.dropout_rate
    opt = config.training.optimizer

    if model_type == "cnn":
        c = m.cnn
        ks = f"k{c.kernel_size[0]}x{c.kernel_size[1]}"
        st = f"s{c.stride[0]}x{c.stride[1]}" if isinstance(c.stride, list) else f"s{c.stride}"
        pad = f"p{c.padding[0]}x{c.padding[1]}" if isinstance(c.padding, list) else f"p{c.padding}"
        pool = f"pool{c.pool_size[0]}x{c.pool_size[1]}_" if c.pool_size else ""
        return f"cnn_{outputs}_{hap}_{tensor_layout}_{wcs}_{norm}_{ks}_f{c.num_filters}_{st}_{pad}_{pool}{hl}_{act}_{dr}_{opt}"

    elif model_type == "cnn2":
        c2 = m.cnn2
        ks1 = c2.kernel_stage1
        return (
            f"cnn2_{outputs}_{hap}_{tensor_layout}_{wcs}_{norm}_"
            f"s1k{ks1[0]}x{ks1[1]}f{c2.num_filters_stage1}_"
            f"s2f{c2.num_filters_stage2}_s3f{c2.num_filters_stage3}_"
            f"gp{c2.global_pool_type}_fc{c2.fc_hidden_size}_"
            f"{hl}_{act}_{dr}_{opt}"
        )

    elif model_type in ("svm", "rf", "xgboost"):
        sk = m.sklearn
        pca = f"pca{sk.pca_components}" if sk.pca_components is not None else "pca_auto"
        if model_type == "svm":
            c_str = str(sk.svm.C).replace(".", "p")
            cal = "cal1" if sk.svm.calibrate_probabilities else "cal0"
            tag = f"svm_C{c_str}_{cal}"
        elif model_type == "rf":
            md = f"md{sk.random_forest.max_depth}" if sk.random_forest.max_depth else "mdNone"
            tag = f"rf_nt{sk.random_forest.n_estimators}_{md}"
        else:
            xgb = sk.xgboost
            lr = str(xgb.learning_rate).replace(".", "p")
            tag = f"xgb_nt{xgb.n_estimators}_md{xgb.max_depth}_lr{lr}"
        return f"{model_type}_{outputs}_{hap}_{tensor_layout}_{wcs}_{norm}_{pca}_{tag}"

    else:
        # NN (default)
        return f"nn_{outputs}_{hap}_{wcs}_{norm}{indel_suffix}_{hl}_{act}_{dr}_{opt}"
