from __future__ import annotations

import torch


def apply_rope(x: torch.Tensor, positions: torch.Tensor, base: float = 10000.0) -> torch.Tensor:
    """Apply RoPE to a tensor shaped (B, H, T, D) using integer genomic positions."""
    if x.shape[-1] % 2 != 0:
        raise ValueError("RoPE requer dimensao par")
    device = x.device
    dtype = x.dtype
    positions = positions.to(device=device, dtype=torch.float32)
    half = x.shape[-1] // 2
    freq_idx = torch.arange(half, device=device, dtype=torch.float32)
    inv_freq = 1.0 / (base ** (freq_idx / half))
    angles = positions[:, None, :, None] * inv_freq[None, None, None, :]
    cos = angles.cos().to(dtype=dtype)
    sin = angles.sin().to(dtype=dtype)
    x_even = x[..., 0::2]
    x_odd = x[..., 1::2]
    out = torch.empty_like(x)
    out[..., 0::2] = x_even * cos - x_odd * sin
    out[..., 1::2] = x_even * sin + x_odd * cos
    return out
