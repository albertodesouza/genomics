from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import torch
from torch import nn

from .model import CNN2Classifier, CNNClassifier


@dataclass
class AttributionResult:
    attribution: torch.Tensor
    target_index: int


class GradCAM:
    def __init__(self, model: nn.Module):
        self.model = model
        self.activations = None
        self.gradients = None
        self._hooks: List = []
        self._register_hooks()

    def _register_hooks(self) -> None:
        if isinstance(self.model, CNNClassifier):
            target_layer = self.model.conv
        elif isinstance(self.model, CNN2Classifier):
            target_layer = self.model.features[4]
        else:
            raise ValueError("GradCAM is supported only for CNN and CNN2 models")

        self._hooks.append(target_layer.register_forward_hook(self._forward_hook))
        self._hooks.append(target_layer.register_full_backward_hook(self._backward_hook))

    def _forward_hook(self, module, inputs, output) -> None:
        self.activations = output.detach()

    def _backward_hook(self, module, grad_inputs, grad_outputs) -> None:
        self.gradients = grad_outputs[0].detach()

    def generate(self, input_tensor: torch.Tensor, target_index: Optional[int] = None) -> AttributionResult:
        was_training = self.model.training
        self.model.eval()

        batch = input_tensor.clone().requires_grad_(True)
        output = self.model(batch)
        if target_index is None:
            target_index = int(output.argmax(dim=1).item())

        self.model.zero_grad()
        one_hot = torch.zeros_like(output)
        one_hot[0, target_index] = 1
        output.backward(gradient=one_hot)

        weights = self.gradients.mean(dim=(2, 3), keepdim=True)
        cam = torch.relu((weights * self.activations).sum(dim=1, keepdim=True))
        cam = torch.nn.functional.interpolate(
            cam,
            size=(input_tensor.shape[1], input_tensor.shape[2]),
            mode="bilinear",
            align_corners=False,
        ).squeeze(0).squeeze(0)

        if was_training:
            self.model.train()
        return AttributionResult(attribution=cam.cpu(), target_index=target_index)

    def close(self) -> None:
        for hook in self._hooks:
            hook.remove()
        self._hooks.clear()


class DeepLIFT:
    def __init__(self, model: nn.Module):
        self.model = model
        self._mean_baseline_cache: Optional[torch.Tensor] = None

    def _get_baseline(self, input_tensor: torch.Tensor, baseline_type: str, dataset=None) -> torch.Tensor:
        if baseline_type == "zeros":
            return torch.zeros_like(input_tensor)
        if baseline_type == "mean":
            if self._mean_baseline_cache is not None:
                return self._mean_baseline_cache.to(input_tensor.device)
            if dataset is None:
                return torch.zeros_like(input_tensor)
            features = []
            for idx in range(min(len(dataset), 1000)):
                feature_tensor, _ = dataset[idx]
                features.append(feature_tensor)
            if not features:
                return torch.zeros_like(input_tensor)
            mean_tensor = torch.stack(features).mean(dim=0, keepdim=True)
            self._mean_baseline_cache = mean_tensor.cpu()
            return mean_tensor.to(input_tensor.device)
        raise ValueError(f"Unsupported baseline_type: {baseline_type}")

    def generate(
        self,
        input_tensor: torch.Tensor,
        target_index: Optional[int] = None,
        baseline_type: str = "zeros",
        dataset=None,
    ) -> AttributionResult:
        was_training = self.model.training
        self.model.eval()

        batch = input_tensor.clone().requires_grad_(True)
        baseline = self._get_baseline(batch, baseline_type=baseline_type, dataset=dataset)
        delta = batch - baseline
        output = self.model(batch)
        if target_index is None:
            target_index = int(output.argmax(dim=1).item())

        self.model.zero_grad()
        one_hot = torch.zeros_like(output)
        one_hot[0, target_index] = 1
        output.backward(gradient=one_hot)
        attribution = (batch.grad.detach() * delta).squeeze(0)

        if was_training:
            self.model.train()
        return AttributionResult(attribution=attribution.cpu(), target_index=target_index)

    def generate_class_mean(
        self,
        target_index: int,
        dataset,
        baseline_type: str = "zeros",
        max_samples: Optional[int] = None,
    ) -> Tuple[Optional[torch.Tensor], Optional[torch.Tensor], int]:
        attributions = []
        inputs = []
        for idx in range(len(dataset)):
            feature_tensor, target_tensor = dataset[idx]
            if int(target_tensor.item()) != target_index:
                continue
            result = self.generate(
                input_tensor=feature_tensor.unsqueeze(0).to(next(self.model.parameters()).device),
                target_index=target_index,
                baseline_type=baseline_type,
                dataset=dataset,
            )
            attributions.append(result.attribution)
            inputs.append(feature_tensor)
            if max_samples is not None and len(attributions) >= max_samples:
                break

        if not attributions:
            return None, None, 0
        return torch.stack(attributions).mean(dim=0), torch.stack(inputs).mean(dim=0), len(attributions)
