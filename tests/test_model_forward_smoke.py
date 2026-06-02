import torch

from genotype_based_predictor.config import PipelineConfig as GenotypeConfig
from genotype_based_predictor.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor
from variant_transformer_predictor.config import ModelConfig
from variant_transformer_predictor.model import VariantTransformerClassifier


def _genotype_config(model_type: str) -> GenotypeConfig:
    return GenotypeConfig.model_validate(
        {
            "dataset_input": {
                "dataset_dir": "/tmp/nonexistent_dataset_for_model_smoke",
                "alphagenome_outputs": ["rna_seq"],
                "haplotype_mode": "H1",
                "tensor_layout": "raw_center_crop",
                "feature_mode": "signals_only",
                "genes_to_use": ["GENE1", "GENE2"],
                "window_center_size": 64,
                "cache_processed_tensors": False,
            },
            "output": {
                "prediction_target": "superpopulation",
                "known_classes": ["A", "B", "C"],
            },
            "model": {
                "type": model_type,
                "hidden_layers": [8],
                "dropout_rate": 0.0,
                "cnn": {
                    "kernel_size": [2, 8],
                    "num_filters": 4,
                    "stride": [1, 4],
                    "padding": 0,
                    "pool_size": None,
                },
                "cnn2": {
                    "num_filters_stage1": 4,
                    "kernel_stage1": [6, 8],
                    "stride_stage1": [6, 4],
                    "num_filters_stage2": 4,
                    "num_filters_stage3": 4,
                    "kernel_stages23": [1, 3],
                    "stride_stages23": [1, 2],
                    "padding_stages23": [0, 1],
                    "fc_hidden_size": 8,
                    "global_pool_type": "avg",
                },
            },
        }
    )


def test_genotype_models_forward_raw_center_crop():
    input_shape = (12, 64)
    x = torch.randn(2, *input_shape)
    for model_cls, model_type in [
        (NNAncestryPredictor, "NN"),
        (CNNAncestryPredictor, "CNN"),
        (CNN2AncestryPredictor, "CNN2"),
    ]:
        model = model_cls(_genotype_config(model_type), input_shape, num_classes=3)
        out = model(x)
        assert out.shape == (2, 3)
        assert torch.isfinite(out).all()


def test_genotype_models_forward_haplotype_channels_input_shape():
    input_shape = (12, 64)
    x = torch.randn(2, 2, 6, 64)
    model = NNAncestryPredictor(_genotype_config("NN"), input_shape, num_classes=3)
    out = model(x)
    assert out.shape == (2, 3)


def test_variant_transformer_forward_smoke():
    cfg = ModelConfig(
        d_type=4,
        d_hap=4,
        d_gene=4,
        d_len=4,
        d_base=4,
        d_allele=8,
        d_model=32,
        layers=1,
        heads=4,
        mlp_ratio=2,
        dropout=0.0,
    )
    model = VariantTransformerClassifier(cfg, num_genes=3, num_classes=5, l_max=4)
    batch = {
        "variant_type": torch.tensor([[0, 1, 2], [1, 0, 2]], dtype=torch.long),
        "haplotype": torch.tensor([[0, 1, 0], [1, 0, 1]], dtype=torch.long),
        "gene": torch.tensor([[0, 1, 2], [2, 1, 0]], dtype=torch.long),
        "length_norm": torch.zeros(2, 3, 1),
        "ref_allele": torch.zeros(2, 3, 4, dtype=torch.long),
        "alt_allele": torch.ones(2, 3, 4, dtype=torch.long),
        "position_relative": torch.tensor([[0.1, 0.2, 0.3], [0.0, 0.4, 0.8]], dtype=torch.float32),
        "attention_mask": torch.ones(2, 4, dtype=torch.bool),
    }

    out = model(batch)
    assert out.shape == (2, 5)
    assert torch.isfinite(out).all()
