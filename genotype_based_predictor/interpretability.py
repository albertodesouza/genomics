# -*- coding: utf-8 -*-
"""
interpretability.py
===================

Ferramentas de interpretabilidade para modelos genômicos.

Classes
-------
GradCAM   : Gradient-weighted Class Activation Mapping (Selvaraju et al., 2017)
DeepLIFT  : Deep Learning Important Features (Shrikumar et al., 2017) — aproximação Rescale Rule

Funções
-------
extract_dna_sequence : Extrai sequência FASTA de uma região genômica
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import torch
import torch.nn as nn
from rich.console import Console

console = Console()


class GradCAM:
    """
    Grad-CAM para CNNs genômicas (CNNAncestryPredictor / CNN2AncestryPredictor).

    Registra hooks na última camada convolucional e calcula a combinação
    ponderada das ativações pelos gradientes médios por canal.

    Parâmetros
    ----------
    model : nn.Module
        Modelo CNN com atributo `conv` (CNN) ou `stage3` (CNN2).
    target_layer : str
        'auto' — detecta automaticamente.
    """

    def __init__(self, model: nn.Module, target_layer: str = "auto"):
        self.model = model
        self.target_layer = target_layer
        self.activations = None
        self.gradients = None
        self._hooks: list = []
        self._setup_hooks()

    def _setup_hooks(self):
        model_type = type(self.model).__name__
        if model_type == "CNNAncestryPredictor":
            target = self.model.conv
            self._layer_name = "conv"
        elif model_type == "CNN2AncestryPredictor":
            # Stage 3 é o último bloco convolucional: Sequential(Conv2d, BN, ReLU)
            # O Conv2d é o primeiro elemento
            target = self.model.stage3[0]
            self._layer_name = "stage3[0]"
        else:
            raise ValueError(f"Grad-CAM não suportado para {model_type}")

        self._hooks.append(target.register_forward_hook(
            lambda m, i, o: setattr(self, "activations", o.detach())
        ))
        self._hooks.append(target.register_full_backward_hook(
            lambda m, gi, go: setattr(self, "gradients", go[0].detach())
        ))
        console.print(f"[cyan]Grad-CAM: camada alvo = {self._layer_name}[/cyan]")

    def generate(self, input_tensor: torch.Tensor, target_class: Optional[int] = None) -> Tuple[torch.Tensor, int]:
        """
        Gera mapa de ativação Grad-CAM.

        Parameters
        ----------
        input_tensor : torch.Tensor
            Tensor [batch, rows, cols].
        target_class : int, optional
            Classe alvo. Se None, usa a predita.

        Returns
        -------
        Tuple[torch.Tensor, int]
            (cam_resized [rows, cols], target_class_idx)
        """
        was_training = self.model.training
        self.model.eval()

        inp = input_tensor.clone().requires_grad_(True)
        output = self.model(inp)

        if target_class is None:
            target_class = output.argmax(dim=1).item()

        self.model.zero_grad()
        n_cls = output.shape[1]
        grad = torch.ones_like(output) / max(n_cls - 1, 1)
        grad[0, target_class] = 0.0
        output.backward(gradient=grad, retain_graph=True)

        weights = self.gradients.mean(dim=(2, 3), keepdim=True)
        cam = torch.relu((weights * self.activations).sum(dim=1, keepdim=True))
        cam = cam.squeeze()

        target_h = len(self.model.rows_to_use) if hasattr(self.model, "rows_to_use") else input_tensor.shape[1]
        target_w = input_tensor.shape[2]

        cam_resized = nn.functional.interpolate(
            cam.unsqueeze(0).unsqueeze(0),
            size=(target_h, target_w),
            mode="bilinear",
            align_corners=False,
        ).squeeze()

        if was_training:
            self.model.train()
        return cam_resized.cpu(), target_class

    def remove_hooks(self):
        for h in self._hooks:
            h.remove()
        self._hooks.clear()


class DeepLIFT:
    """
    DeepLIFT via aproximação Rescale Rule: atribuição = gradiente × (input − baseline).

    Parâmetros
    ----------
    model : nn.Module
        Qualquer modelo (NN, CNN, CNN2).
    """

    def __init__(self, model: nn.Module):
        self.model = model
        self._baseline_cache: Optional[torch.Tensor] = None
        self._class_mean_cache: Dict[int, torch.Tensor] = {}
        self._class_input_mean_cache: Dict[int, Tuple[torch.Tensor, int]] = {}

    def _get_baseline(self, input_tensor: torch.Tensor, baseline_type: str, dataset=None) -> torch.Tensor:
        if baseline_type == "zeros":
            return torch.zeros_like(input_tensor)
        elif baseline_type == "mean":
            if self._baseline_cache is not None:
                return self._baseline_cache.to(input_tensor.device)
            if dataset is None:
                return torch.zeros_like(input_tensor)
            console.print("[cyan]Calculando baseline (média do dataset)...[/cyan]")
            samples = [dataset[i][0] for i in range(min(len(dataset), 10000))]
            self._baseline_cache = torch.stack(samples).mean(dim=0)
            return self._baseline_cache.unsqueeze(0).to(input_tensor.device)
        raise ValueError(f"baseline_type inválido: {baseline_type}")

    def generate(self, input_tensor: torch.Tensor, target_class: Optional[int] = None,
                 baseline_type: str = "zeros", dataset=None) -> Tuple[torch.Tensor, int]:
        """
        Gera atribuições DeepLIFT.

        Returns
        -------
        Tuple[torch.Tensor, int]
            (atribuições [rows, cols], target_class_idx)
        """
        was_training = self.model.training
        self.model.eval()

        baseline = self._get_baseline(input_tensor, baseline_type, dataset)
        delta = input_tensor - baseline

        inp = input_tensor.clone().requires_grad_(True)
        output = self.model(inp)

        if target_class is None:
            target_class = output.argmax(dim=1).item()

        self.model.zero_grad()
        one_hot = torch.zeros_like(output)
        one_hot[0, target_class] = 1.0
        output.backward(gradient=one_hot)

        attributions = (inp.grad.detach() * delta).squeeze(0)

        if was_training:
            self.model.train()
        return attributions.cpu(), target_class

    def generate_class_mean(self, target_class_idx: int, dataset: Any,
                            baseline_type: str = "zeros") -> Tuple[Optional[torch.Tensor], Optional[torch.Tensor], int]:
        """
        Calcula a média das atribuições DeepLIFT para todas as amostras de uma classe.

        Returns
        -------
        Tuple[mean_attributions, mean_input, num_samples]
        """
        if target_class_idx in self._class_mean_cache:
            mean_input, n = self._class_input_mean_cache[target_class_idx]
            return self._class_mean_cache[target_class_idx], mean_input, n

        device = next(self.model.parameters()).device
        samples = [dataset[i][0] for i in range(len(dataset)) if dataset[i][1] == target_class_idx]
        if not samples:
            return None, None, 0

        mean_input = torch.stack(samples).mean(dim=0)
        attrs = []
        for s in samples:
            a, _ = self.generate(s.unsqueeze(0).to(device), target_class=target_class_idx,
                                 baseline_type=baseline_type, dataset=dataset)
            attrs.append(a)

        mean_attr = torch.stack(attrs).mean(dim=0)
        self._class_mean_cache[target_class_idx] = mean_attr
        self._class_input_mean_cache[target_class_idx] = (mean_input, len(samples))
        return mean_attr, mean_input, len(samples)


def extract_dna_sequence(dataset_dir: Path, sample_id: str, gene_name: str,
                         center_position: int, window_center_size: int,
                         sequence_length: int = 1000, haplotype: str = "H1") -> Optional[Tuple[str, str]]:
    """
    Extrai sequência de DNA centrada numa posição genômica de um arquivo FASTA individual.

    Parameters
    ----------
    dataset_dir : Path
        Diretório raiz do dataset.
    sample_id : str
        ID do indivíduo (ex: HG02635).
    gene_name : str
        Nome do gene (ex: DDB1).
    center_position : int
        Posição genômica central (para o header FASTA).
    window_center_size : int
        Tamanho da janela central usada no processamento.
    sequence_length : int
        Comprimento da sequência a extrair em bp.
    haplotype : str
        'H1' ou 'H2'.

    Returns
    -------
    Optional[Tuple[str, str]]
        (header, sequence) ou None se não encontrado.
    """
    fasta_path = (Path(dataset_dir) / "individuals" / sample_id / "windows" / gene_name
                  / f"{sample_id}.{haplotype}.window.fixed.fa")
    if not fasta_path.exists():
        console.print(f"[yellow]⚠ FASTA não encontrado: {fasta_path}[/yellow]")
        return None
    try:
        with open(fasta_path) as f:
            lines = f.readlines()
        if len(lines) < 2:
            return None
        sequence = "".join(l.strip() for l in lines[1:])
        center = len(sequence) // 2
        half = sequence_length // 2
        seq = sequence[max(0, center - half): min(len(sequence), center + half)]
        header = f">{sample_id}_{haplotype}_{gene_name}_center_{center_position}"
        return header, seq
    except Exception as e:
        console.print(f"[yellow]⚠ Erro ao ler FASTA: {e}[/yellow]")
        return None
