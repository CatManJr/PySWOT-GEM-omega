"""Gaussian smoothing matching MATLAB ``smoothdata(..., 'gaussian', window)``."""

from __future__ import annotations

import numpy as np
from numba import njit

from swot_qg_omega.constants import GAUSS_WINDOW


def _gaussian_window_weights(window: int) -> tuple[np.ndarray, int, int]:
    """
    MATLAB ``smoothdata`` gaussian: ``sigma = window/5``, truncated to ``window`` taps.

    For even ``window``, the window is centered on the current and previous
    elements (``back = window//2``, ``forward = window//2 - 1``).
    """
    sigma = window / 5.0
    back = window // 2
    forward = window - back - 1
    offsets = np.arange(-back, forward + 1, dtype=np.float64)
    weights = np.exp(-0.5 * (offsets / sigma) ** 2)
    return weights, back, forward


@njit(cache=True)
def _smooth_lines_numba(
    flat: np.ndarray,
    weights: np.ndarray,
    back: int,
    forward: int,
) -> np.ndarray:
    n_rows, n = flat.shape
    out = np.empty_like(flat)
    for row in range(n_rows):
        line = flat[row]
        for i in range(n):
            i0 = i - back
            if i0 < 0:
                i0 = 0
            i1 = i + forward + 1
            if i1 > n:
                i1 = n
            w0 = back - (i - i0)
            num = 0.0
            den = 0.0
            for j in range(i1 - i0):
                val = line[i0 + j]
                if np.isfinite(val):
                    ww = weights[w0 + j]
                    num += ww * val
                    den += ww
            if den > 0.0:
                out[row, i] = num / den
            else:
                out[row, i] = np.nan
    return out


def _smoothdata_gaussian_1d(
    field: np.ndarray,
    axis: int,
    window: int = GAUSS_WINDOW,
) -> np.ndarray:
    """Apply MATLAB-like ``smoothdata(..., dim, 'gaussian', window)`` along one axis."""
    weights, back, forward = _gaussian_window_weights(window)
    moved = np.moveaxis(np.asarray(field, dtype=np.float64), axis, -1)
    flat = np.ascontiguousarray(moved.reshape(-1, moved.shape[-1]))
    out_flat = _smooth_lines_numba(flat, weights, back, forward)
    return np.moveaxis(out_flat.reshape(moved.shape), -1, axis)


def smooth_gaussian_2d(field: np.ndarray, window: int = GAUSS_WINDOW) -> np.ndarray:
    """
    Match CreateDensityFields.m:

    ``smoothdata(X, 1, 'gaussian', window);`` then ``smoothdata(..., 2, 'gaussian', window)``.
    """
    out = _smoothdata_gaussian_1d(field, axis=0, window=window)
    out = _smoothdata_gaussian_1d(out, axis=1, window=window)
    return out


def smooth_depth_stack(stack: np.ndarray, window: int = GAUSS_WINDOW) -> np.ndarray:
    """Smooth each depth level of an array with trailing depth axis."""
    out = np.empty_like(stack, dtype=float)
    for k in range(stack.shape[-1]):
        out[..., k] = smooth_gaussian_2d(stack[..., k], window=window)
    return out
