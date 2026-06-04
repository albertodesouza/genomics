from __future__ import annotations

import argparse
from pathlib import Path

import torch
from torch.utils.data import DataLoader

from .config import get_experiment_dir, load_config
from .dataset import VariantTokenDataset, collate_variant_tokens
from .evaluation import evaluate, save_results
from .model import VariantTransformerClassifier


def main() -> int:
    parser = argparse.ArgumentParser(description="Avalia checkpoint do Variant Transformer")
    parser.add_argument("config_path")
    parser.add_argument("--checkpoint", default="best_accuracy")
    parser.add_argument("--split", choices=["train", "val", "test"], default="test")
    parser.add_argument("--experiment-dir", type=Path, default=None)
    args = parser.parse_args()
    config = load_config(Path(args.config_path))
    experiment_dir = args.experiment_dir or get_experiment_dir(config)
    ds = VariantTokenDataset(
        config.dataset.processed_dir,
        split=args.split,
        max_sequence_length=config.dataset.max_sequence_length,
        truncate_policy=config.dataset.truncate_policy,
        loading_strategy=config.dataset.loading_strategy,
    )
    loader = DataLoader(ds, batch_size=config.training.batch_size, shuffle=False, num_workers=config.training.num_workers, collate_fn=collate_variant_tokens)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = VariantTransformerClassifier(config.model, int(ds.metadata["num_genes"]), ds.get_num_classes(), int(ds.metadata["l_max"])).to(device)
    checkpoint = Path(args.checkpoint)
    if not checkpoint.exists():
        checkpoint = experiment_dir / "models" / (args.checkpoint if args.checkpoint.endswith(".pt") else f"{args.checkpoint}.pt")
    state = torch.load(checkpoint, map_location=device, weights_only=True)
    model.load_state_dict(state.get("model_state_dict", state))
    results = evaluate(model, loader, ds.get_class_names(), device, f"{args.split}_{checkpoint.stem}")
    if results:
        save_results(results, experiment_dir / f"{args.split}_{checkpoint.stem}_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
