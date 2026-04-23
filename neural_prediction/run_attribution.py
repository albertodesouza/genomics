#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from .config import load_config
from .data import create_dataloaders
from .interpretability import DeepLIFT, GradCAM
from .model import build_model
from .train import load_checkpoint_if_available


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run attribution methods for neural_prediction models")
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--method", choices=["gradcam", "deeplift"], required=True)
    parser.add_argument("--sample-index", type=int, default=0)
    parser.add_argument("--output", required=True)
    parser.add_argument("--baseline", default="zeros", choices=["zeros", "mean"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(Path(args.config))
    data_bundle = create_dataloaders(config)
    model = build_model(config, data_bundle.dataset.get_input_shape(), data_bundle.dataset.get_num_classes())
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    load_checkpoint_if_available(model, args.checkpoint, device)

    feature_tensor, target_tensor = data_bundle.dataset[args.sample_index]
    input_tensor = feature_tensor.unsqueeze(0).to(device)

    if args.method == "gradcam":
        runner = GradCAM(model)
        result = runner.generate(input_tensor)
        runner.close()
    else:
        runner = DeepLIFT(model)
        result = runner.generate(input_tensor, baseline_type=args.baseline, dataset=data_bundle.dataset)

    payload = {
        "sample_index": args.sample_index,
        "target_index": int(result.target_index),
        "target_value": int(target_tensor.item()) if target_tensor.ndim == 0 else None,
        "shape": list(result.attribution.shape),
        "attribution": result.attribution.tolist(),
    }
    with Path(args.output).open("w", encoding="utf-8") as handle:
        json.dump(payload, handle)
        handle.write("\n")


if __name__ == "__main__":
    main()
