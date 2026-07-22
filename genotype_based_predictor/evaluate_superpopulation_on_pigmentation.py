from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import torch
from rich.console import Console
from rich.table import Table
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_fscore_support

from genotype_based_predictor.config import generate_experiment_name, load_config
from genotype_based_predictor.data_pipeline import prepare_data
from genotype_based_predictor.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor

console = Console()


def _build_model(config, dataset, num_classes: int):
    input_shape = dataset.get_input_shape()
    model_type = config.model.type.upper()
    if model_type == "NN":
        return NNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN":
        return CNNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN2":
        return CNN2AncestryPredictor(config, input_shape, num_classes)
    raise ValueError(f"Modelo PyTorch nao suportado: {config.model.type}")


def _resolve_checkpoint(experiment_dir: Path, checkpoint: str) -> Path:
    path = Path(checkpoint)
    if path.exists():
        return path.resolve()
    if path.suffix != ".pt":
        path = path.with_suffix(".pt")
    return experiment_dir / "models" / path.name


def _select_loader(split_name: str, train_loader, val_loader, test_loader):
    split = split_name.lower()
    if split == "train":
        return train_loader
    if split == "val":
        return val_loader
    if split == "test":
        return test_loader
    raise ValueError("--split deve ser train, val ou test")


def _parse_mapping(mapping_items: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for item in mapping_items:
        if "=" not in item:
            raise ValueError(f"Mapeamento invalido: {item}. Use SUPERPOP=PIGMENTATION")
        left, right = item.split("=", 1)
        mapping[left.strip()] = right.strip()
    return mapping


def _get_sample_id(dataset, idx: int) -> str:
    if hasattr(dataset, "get_sample_id"):
        return dataset.get_sample_id(idx)
    if hasattr(dataset, "dataset") and hasattr(dataset, "indices"):
        return _get_sample_id(dataset.dataset, dataset.indices[idx])
    raise TypeError(f"Dataset sem get_sample_id: {type(dataset).__name__}")


def _collect_targets_by_sample_id(loader) -> Dict[str, int]:
    dataset = loader.dataset
    targets_by_sample_id: Dict[str, int] = {}
    for idx in range(len(dataset)):
        item = dataset[idx]
        _features, target = item[:2]
        sample_id = _get_sample_id(dataset, idx)
        targets_by_sample_id[sample_id] = int(target.item())
    return targets_by_sample_id


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate a superpopulation checkpoint on a binary pigmentation split by mapping predicted superpopulations."
    )
    parser.add_argument("superpopulation_config", type=Path, help="YAML usado para treinar o modelo de superpopulacao")
    parser.add_argument("pigmentation_config", type=Path, help="YAML que define o split/target de pigmentacao")
    parser.add_argument("--checkpoint", default="best_accuracy", help="Checkpoint alias/path do experimento de superpopulacao")
    parser.add_argument("--split", default="test", choices=["train", "val", "test"], help="Split da config de pigmentacao")
    parser.add_argument("--superpopulation-experiment-dir", type=Path, default=None, help="Diretorio do experimento de superpopulacao")
    parser.add_argument("--pigmentation-experiment-dir", type=Path, default=None, help="Diretorio usado para carregar/criar o cache de pigmentacao")
    parser.add_argument("--output-name", default=None, help="Nome do JSON de saida sem .json")
    parser.add_argument(
        "--map",
        action="append",
        default=["AFR=strong pigmentation", "EUR=weak pigmentation"],
        help="Mapeamento SUPERPOP=PIGMENTATION. Pode repetir. Default: AFR/EUR",
    )
    args = parser.parse_args()

    super_config = load_config(args.superpopulation_config.resolve())
    pigmentation_config = load_config(args.pigmentation_config.resolve())

    super_experiment_dir = args.superpopulation_experiment_dir or (
        Path(super_config.dataset_input.processed_cache_dir) / generate_experiment_name(super_config)
    )
    pigmentation_experiment_dir = args.pigmentation_experiment_dir or (
        Path(pigmentation_config.dataset_input.processed_cache_dir) / generate_experiment_name(pigmentation_config)
    )
    super_experiment_dir = super_experiment_dir.resolve()
    pigmentation_experiment_dir = pigmentation_experiment_dir.resolve()
    checkpoint_path = _resolve_checkpoint(super_experiment_dir, args.checkpoint)
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint nao encontrado: {checkpoint_path}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    console.print(f"[green]Experimento superpopulacao:[/green] {super_experiment_dir}")
    console.print(f"[green]Experimento/split pigmentacao:[/green] {pigmentation_experiment_dir}")
    console.print(f"[green]Checkpoint:[/green] {checkpoint_path}")
    console.print(f"[green]Split:[/green] {args.split}")
    console.print(f"[green]Device:[/green] {device}")

    pigmentation_full_ds, pigmentation_train_loader, pigmentation_val_loader, pigmentation_test_loader = prepare_data(
        pigmentation_config, pigmentation_experiment_dir
    )
    pigmentation_loader = _select_loader(args.split, pigmentation_train_loader, pigmentation_val_loader, pigmentation_test_loader)
    pigmentation_targets_by_sample_id = _collect_targets_by_sample_id(pigmentation_loader)

    super_full_ds, super_train_loader, super_val_loader, super_test_loader = prepare_data(super_config, super_experiment_dir)

    super_classes = super_config.output.known_classes or ["AFR", "AMR", "EAS", "EUR", "SAS"]
    model = _build_model(super_config, super_full_ds, len(super_classes)).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint.get("model_state_dict", checkpoint))
    model.eval()

    super_idx_to_target = {idx: target for idx, target in enumerate(super_classes)}
    super_to_pigmentation = _parse_mapping(args.map)
    pigmentation_target_to_idx = pigmentation_full_ds.target_to_idx
    pigmentation_idx_to_target = pigmentation_full_ds.idx_to_target

    all_preds: List[int] = []
    all_targets: List[int] = []
    skipped_unmapped_predictions = 0
    seen_sample_ids = set()
    source_split_counts = {"train": 0, "val": 0, "test": 0}

    with torch.no_grad():
        for source_split, loader in [
            ("train", super_train_loader),
            ("val", super_val_loader),
            ("test", super_test_loader),
        ]:
            dataset = loader.dataset
            for idx in range(len(dataset)):
                sample_id = _get_sample_id(dataset, idx)
                target_idx = pigmentation_targets_by_sample_id.get(sample_id)
                if target_idx is None:
                    continue

                item = dataset[idx]
                features = item[0].unsqueeze(0).to(device, non_blocking=True)
                outputs = model(features)
                pred_super = super_idx_to_target[int(outputs.argmax(dim=1).item())]
                pred_pigmentation = super_to_pigmentation.get(pred_super)
                seen_sample_ids.add(sample_id)
                source_split_counts[source_split] += 1
                if pred_pigmentation is None:
                    skipped_unmapped_predictions += 1
                    continue
                if pred_pigmentation not in pigmentation_target_to_idx:
                    raise ValueError(f"Classe de pigmentacao mapeada nao existe no dataset: {pred_pigmentation}")
                all_preds.append(pigmentation_target_to_idx[pred_pigmentation])
                all_targets.append(int(target_idx))

    missing_in_superpopulation_cache = len(set(pigmentation_targets_by_sample_id) - seen_sample_ids)

    if not all_targets:
        raise RuntimeError("Nenhuma amostra valida sobrou apos aplicar o mapeamento de superpopulacao para pigmentacao")

    labels = sorted(pigmentation_idx_to_target)
    target_names = [pigmentation_idx_to_target[i] for i in labels]
    precision, recall, f1, _ = precision_recall_fscore_support(
        all_targets, all_preds, average="weighted", zero_division=0
    )
    accuracy = accuracy_score(all_targets, all_preds)
    cm = confusion_matrix(all_targets, all_preds, labels=labels)
    report = classification_report(all_targets, all_preds, labels=labels, target_names=target_names, zero_division=0)

    results: Dict[str, Any] = {
        "accuracy": float(accuracy),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "confusion_matrix": cm.tolist(),
        "classification_report": report,
        "num_samples": len(all_targets),
        "skipped_unmapped_predictions": skipped_unmapped_predictions,
        "missing_in_superpopulation_cache": missing_in_superpopulation_cache,
        "pigmentation_split_samples": len(pigmentation_targets_by_sample_id),
        "superpopulation_source_split_counts": source_split_counts,
        "superpopulation_classes": super_classes,
        "superpopulation_to_pigmentation": super_to_pigmentation,
    }

    table = Table(title="Superpopulation checkpoint on pigmentation split", show_header=True)
    table.add_column("Metrica")
    table.add_column("Valor", justify="right")
    table.add_row("Accuracy", f"{accuracy:.4f}")
    table.add_row("Precision (weighted)", f"{precision:.4f}")
    table.add_row("Recall (weighted)", f"{recall:.4f}")
    table.add_row("F1 (weighted)", f"{f1:.4f}")
    table.add_row("Amostras avaliadas", str(len(all_targets)))
    table.add_row("Predicoes descartadas", str(skipped_unmapped_predictions))
    table.add_row("Ausentes no cache superpop", str(missing_in_superpopulation_cache))
    table.add_row("Origem no treino superpop", json.dumps(source_split_counts, sort_keys=True))
    console.print(table)
    console.print(report)

    output_name = args.output_name or f"{args.split}_{checkpoint_path.stem}_superpopulation_mapped_to_pigmentation"
    output_path = pigmentation_experiment_dir / f"{output_name}.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    console.print(f"[green]Resultados salvos:[/green] {output_path}")


if __name__ == "__main__":
    main()
