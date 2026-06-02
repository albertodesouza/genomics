from __future__ import annotations

from typing import Tuple

import numpy as np


def block_reduce_2d(arr: np.ndarray, target_shape: Tuple[int, int], func: str = "max") -> np.ndarray:
    """Resize a 2D array with interpolation for upscaling and block reduction for downscaling."""
    h, w = arr.shape
    th, tw = target_shape

    agg_map = {"max": np.max, "min": np.min, "mean": np.mean}
    if func not in agg_map:
        raise ValueError(f"Função de agregação desconhecida: '{func}'. Use 'max', 'min' ou 'mean'.")
    agg_func = agg_map[func]

    if tw >= w:
        from scipy import ndimage
        temp = ndimage.zoom(arr, (1, tw / w), order=1)
    else:
        bw = w / tw
        temp = np.zeros((h, tw), dtype=arr.dtype)
        for j in range(tw):
            x_start = int(j * bw)
            x_end = min(int((j + 1) * bw), w)
            if x_end > x_start:
                temp[:, j] = agg_func(arr[:, x_start:x_end], axis=1)

    h_temp = temp.shape[0]
    if th >= h_temp:
        from scipy import ndimage
        return ndimage.zoom(temp, (th / h_temp, 1), order=1)

    bh = h_temp / th
    result = np.zeros((th, tw), dtype=arr.dtype)
    for i in range(th):
        y_start = int(i * bh)
        y_end = min(int((i + 1) * bh), h_temp)
        if y_end > y_start:
            result[i, :] = agg_func(temp[y_start:y_end, :], axis=0)
    return result
