"""Overlap scores used to compare this pipeline with the published NetCDF."""

from __future__ import annotations

import numpy as np


def overlap_corr(reference: np.ndarray, estimate: np.ndarray) -> float:
    """Pearson correlation on the finite overlap of two arrays."""
    mask = np.isfinite(reference) & np.isfinite(estimate)
    if int(mask.sum()) < 10:
        return float("nan")
    left = reference[mask]
    right = estimate[mask]
    if left.std() == 0.0 or right.std() == 0.0:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def to_cross_along_depth(values: np.ndarray) -> np.ndarray:
    """
    Published NetCDF stores full-swath fields as (depth, along, cross).

    Return (cross, along, depth), which is the layout written by this package.
    """
    array = np.asarray(values, dtype=float)
    if array.ndim == 3 and array.shape[0] in (245, 490) and array.shape[-1] in (52, 69):
        return np.transpose(array, (2, 1, 0))
    return array
