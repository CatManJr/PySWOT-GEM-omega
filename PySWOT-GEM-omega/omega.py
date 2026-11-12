"""Quasi-geostrophic omega inversion on tiled swath domains via xinvert."""

from __future__ import annotations

import numpy as np
import xarray as xr
from xinvert import invert_omega

from swot_qg_omega.constants import GRID_SPACING_M

# Match the paper workflow: Dirichlet in z, periodic horizontal tile, SOR.
DEFAULT_MX_LOOP = 20000
DEFAULT_TOLERANCE = 1e-7


def invert_tiled_omega(
    qf_trunc: np.ndarray,
    n2_trunc: np.ndarray,
    depth_trunc: np.ndarray,
    f0: float,
    *,
    mx_loop: int = DEFAULT_MX_LOOP,
    tolerance: float = DEFAULT_TOLERANCE,
) -> np.ndarray:
    """
    Invert $F$ with ``xinvert.invert_omega`` (paper / Zenodo workflow).

    Solves
    $\\overline{N^{2}}\\nabla_{h}^{2}w+f^{2}w_{zz}=F$
    on the mirror-tiled Cartesian swath. ``qf_trunc`` has shape
    (cross, along, depth). Returns $w$ in $\\mathrm{m\\,s}^{-1}$.

    Horizontal BCs request periodic–periodic as in the paper tiling; current
    ``xinvert`` SOR kernels only enforce periodic on the last axis (along-track
    here). Top/bottom use Dirichlet $w=0$ via a zero ``icbc``.
    """
    mask = np.isfinite(qf_trunc)
    qf = np.where(mask, qf_trunc, 0.0)
    n2 = np.asarray(n2_trunc, dtype=float)
    n2 = np.where(np.isfinite(n2) & (n2 > 0), n2, np.nanmean(n2[n2 > 0]))

    n_cross, n_along, _ = qf.shape
    cross_m = np.arange(n_cross) * GRID_SPACING_M
    along_m = np.arange(n_along) * GRID_SPACING_M
    # (depth, along, cross): cross-track last. xinvert only enforces periodic
    # BC on the last axis; the mirror tile is narrowest in cross-track.
    F = xr.DataArray(
        np.transpose(qf, (2, 1, 0)),
        dims=["lev", "y", "x"],
        coords={"lev": depth_trunc, "y": along_m, "x": cross_m},
        name="QF",
    )
    N2 = xr.DataArray(n2, dims=["lev"], coords={"lev": depth_trunc})
    icbc = xr.zeros_like(F)

    w = invert_omega(
        F,
        dims=["lev", "y", "x"],
        coords="cartesian",
        icbc=icbc,
        mParams={"f0": float(f0), "beta": 0.0, "N2": N2},
        iParams={
            "BCs": ["fixed", "periodic", "periodic"],
            "mxLoop": mx_loop,
            "tolerance": tolerance,
            "optArg": None,
        },
    )
    # xinvert returns (lev, y=along, x=cross) -> (cross, along, depth)
    w_np = np.transpose(np.asarray(w.values, dtype=float), (2, 1, 0))
    w_np[~mask] = np.nan
    return w_np


def until_swath(
    w_left: np.ndarray,
    w_right: np.ndarray,
    cross_size: int,
    along_size: int,
    orientation: int,
) -> np.ndarray:
    """Place left/right tiled ``xinvert`` solutions back onto the swath grid."""
    n_z = w_left.shape[-1]
    wfull = np.full((cross_size, along_size, n_z), np.nan, dtype=float)
    if orientation == 1:
        wfull[4:31, :, :] = w_left[:27, :along_size, :]
        wfull[39:66, :, :] = w_right[:27, :along_size, :]
    else:
        wfull[4:30, :, :] = w_left[:26, :along_size, :]
        wfull[39:65, :, :] = w_right[:26, :along_size, :]
    return wfull
