from __future__ import annotations

import math
from typing import Tuple

import torch
from torch import nn

from .config import PredictionConfig


class MLPClassifier(nn.Module):
    def __init__(self, input_shape: Tuple[int, int], hidden_layers: list[int], dropout_rate: float, num_classes: int):
        super().__init__()
        input_size = input_shape[0] * input_shape[1]
        layers = [nn.Flatten()]
        current_size = input_size
        for hidden_size in hidden_layers:
            layers.append(nn.Linear(current_size, hidden_size))
            layers.append(nn.ReLU())
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            current_size = hidden_size
        layers.append(nn.Linear(current_size, num_classes))
        self.network = nn.Sequential(*layers)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.network(inputs)


class CNNClassifier(nn.Module):
    def __init__(self, config: PredictionConfig, input_shape: Tuple[int, int], num_classes: int):
        super().__init__()
        num_rows, effective_size = input_shape
        cnn_config = config.model.cnn or {}
        kernel_size = tuple(cnn_config.get("kernel_size", [3, 9]))
        stride = tuple(cnn_config.get("stride", [1, 2]))
        padding = cnn_config.get("padding", 0)
        pool_size = cnn_config.get("pool_size")
        num_filters = int(cnn_config.get("num_filters", 16))

        self.activation = _build_activation(config.model.activation)
        self.conv = nn.Conv2d(1, num_filters, kernel_size=kernel_size, stride=stride, padding=padding)
        self.pool = nn.MaxPool2d(tuple(pool_size)) if pool_size is not None else None

        conv_out_h, conv_out_w = _conv2d_output_shape(num_rows, effective_size, kernel_size, stride, padding)
        if self.pool is not None:
            pool_kernel = tuple(pool_size)
            conv_out_h = conv_out_h // pool_kernel[0]
            conv_out_w = conv_out_w // pool_kernel[1]

        flattened_size = num_filters * conv_out_h * conv_out_w
        layers = []
        prev_size = flattened_size
        for hidden_size in config.model.hidden_layers:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(_build_activation(config.model.activation))
            if config.model.dropout_rate > 0:
                layers.append(nn.Dropout(config.model.dropout_rate))
            prev_size = hidden_size
        layers.append(nn.Linear(prev_size, num_classes))
        self.classifier = nn.Sequential(*layers)
        self._initialize_weights(config.model.activation)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        x = inputs.unsqueeze(1)
        x = self.activation(self.conv(x))
        if self.pool is not None:
            x = self.pool(x)
        x = x.reshape(x.size(0), -1)
        return self.classifier(x)

    def _initialize_weights(self, activation: str) -> None:
        for module in self.modules():
            if isinstance(module, nn.Conv2d):
                _init_conv_or_linear(module, activation, is_output=False)
            elif isinstance(module, nn.Linear):
                is_output = module.out_features == self.classifier[-1].out_features
                _init_conv_or_linear(module, activation, is_output=is_output)


class CNN2Classifier(nn.Module):
    def __init__(self, config: PredictionConfig, input_shape: Tuple[int, int], num_classes: int):
        super().__init__()
        num_rows, effective_size = input_shape
        cnn2_config = config.model.cnn2 or {}
        stage1_filters = int(cnn2_config.get("num_filters_stage1", 16))
        stage2_filters = int(cnn2_config.get("num_filters_stage2", 32))
        stage3_filters = int(cnn2_config.get("num_filters_stage3", 64))
        kernel_stage1 = tuple(cnn2_config.get("kernel_stage1", [min(num_rows, 6), 32]))
        stride_stage1 = tuple(cnn2_config.get("stride_stage1", [max(1, kernel_stage1[0]), 16]))
        kernel_stages23 = tuple(cnn2_config.get("kernel_stages23", [1, 5]))
        stride_stages23 = tuple(cnn2_config.get("stride_stages23", [1, 2]))
        padding_stages23 = tuple(cnn2_config.get("padding_stages23", [0, 2]))
        fc_hidden_size = int(cnn2_config.get("fc_hidden_size", 128))
        pool_type = str(cnn2_config.get("global_pool_type", "max")).lower()

        self.features = nn.Sequential(
            nn.Conv2d(1, stage1_filters, kernel_size=kernel_stage1, stride=stride_stage1),
            nn.ReLU(),
            nn.Conv2d(stage1_filters, stage2_filters, kernel_size=kernel_stages23, stride=stride_stages23, padding=padding_stages23),
            nn.ReLU(),
            nn.Conv2d(stage2_filters, stage3_filters, kernel_size=kernel_stages23, stride=stride_stages23, padding=padding_stages23),
            nn.ReLU(),
        )

        conv1_h, conv1_w = _conv2d_output_shape(num_rows, effective_size, kernel_stage1, stride_stage1, 0)
        conv2_h, conv2_w = _conv2d_output_shape(conv1_h, conv1_w, kernel_stages23, stride_stages23, padding_stages23)
        conv3_h, conv3_w = _conv2d_output_shape(conv2_h, conv2_w, kernel_stages23, stride_stages23, padding_stages23)
        pool_kernel = (1, max(1, conv3_w))
        self.global_pool = nn.MaxPool2d(pool_kernel) if pool_type == "max" else nn.AvgPool2d(pool_kernel)

        flattened_size = stage3_filters * max(1, conv3_h)
        self.classifier = nn.Sequential(
            nn.Linear(flattened_size, fc_hidden_size),
            nn.ReLU(),
            nn.Dropout(config.model.dropout_rate),
            nn.Linear(fc_hidden_size, num_classes),
        )
        self._initialize_weights()

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        x = inputs.unsqueeze(1)
        x = self.features(x)
        x = self.global_pool(x)
        x = x.reshape(x.size(0), -1)
        return self.classifier(x)

    def _initialize_weights(self) -> None:
        for module in self.modules():
            if isinstance(module, nn.Conv2d):
                nn.init.kaiming_normal_(module.weight, mode="fan_in", nonlinearity="relu")
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Linear):
                is_output = module.out_features == self.classifier[-1].out_features
                _init_conv_or_linear(module, "relu", is_output=is_output)


def _build_activation(name: str) -> nn.Module:
    lowered = name.lower()
    if lowered == "relu":
        return nn.ReLU()
    if lowered == "tanh":
        return nn.Tanh()
    if lowered == "sigmoid":
        return nn.Sigmoid()
    raise ValueError(f"Unsupported activation: {name}")


def _conv2d_output_shape(height: int, width: int, kernel_size: Tuple[int, int], stride: Tuple[int, int], padding) -> Tuple[int, int]:
    if isinstance(padding, int):
        pad_h = pad_w = padding
    else:
        pad_h, pad_w = padding
    out_h = math.floor((height + 2 * pad_h - kernel_size[0]) / stride[0] + 1)
    out_w = math.floor((width + 2 * pad_w - kernel_size[1]) / stride[1] + 1)
    if out_h <= 0 or out_w <= 0:
        raise ValueError(
            "Invalid convolution configuration for input shape "
            f"({height}, {width}) with kernel={kernel_size}, stride={stride}, padding={padding}"
        )
    return out_h, out_w


def _init_conv_or_linear(module: nn.Module, activation: str, is_output: bool) -> None:
    if is_output:
        nn.init.xavier_normal_(module.weight)
    elif activation.lower() == "relu":
        nn.init.kaiming_normal_(module.weight, mode="fan_in", nonlinearity="relu")
    else:
        nn.init.xavier_normal_(module.weight)
    if module.bias is not None:
        nn.init.zeros_(module.bias)


def build_model(config: PredictionConfig, input_shape: Tuple[int, int], num_classes: int) -> nn.Module:
    model_type = config.model.type.upper()
    if model_type in {"MLP", "NN"}:
        return MLPClassifier(
            input_shape=input_shape,
            hidden_layers=config.model.hidden_layers,
            dropout_rate=float(config.model.dropout_rate),
            num_classes=num_classes,
        )
    if model_type == "CNN":
        return CNNClassifier(config=config, input_shape=input_shape, num_classes=num_classes)
    if model_type == "CNN2":
        return CNN2Classifier(config=config, input_shape=input_shape, num_classes=num_classes)

        raise ValueError(f"Unsupported model.type for neural_prediction: {model_type}")
