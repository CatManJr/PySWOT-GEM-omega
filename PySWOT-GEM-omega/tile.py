"""Doubly periodic tiling of left and right KaRIn swaths for elliptic inversion."""

from __future__ import annotations

import numpy as np

from swot_qg_omega.constants import N_DEPTH, SWATH_EDGES_ORIENT1, SWATH_EDGES_ORIENT2


def swath_edges(orientation: int) -> dict[str, tuple[int, int]]:
    """Return half-open Python slices (i0, i1) for left and right swaths."""
    src = SWATH_EDGES_ORIENT1 if orientation == 1 else SWATH_EDGES_ORIENT2
    return {
        "left": (src["left"][0] - 1, src["left"][1]),
        "right": (src["right"][0] - 1, src["right"][1]),
    }


def seal_periodic_cross_seam(tile: np.ndarray) -> np.ndarray:
    """
    Match CreateDensityFields.m cleanup used for pass 155 (and most other swaths):

    ``tile(1,:,:) = (tile(2,:,:) + tile(end,:,:)) / 2``
    so the cross-track periodic seam is continuous before omega inversion.
    """
    out = np.array(tile, dtype=float, copy=True)
    out[0, :, :] = 0.5 * (out[1, :, :] + out[-1, :, :])
    return out


def mirror_tile_horizontal(block: np.ndarray, *, seal_seam: bool = True) -> np.ndarray:
    """
    Reflect a swath block in cross-track and along-track to form a periodic tile.

    Input shape is (n_cross, n_along, n_z).
    """
    n_cross, n_along, n_z = block.shape
    tilesize_lon = (n_cross - 1) * 2
    tilesize_lat = 2 * n_along - 2
    tile = np.full((tilesize_lon, tilesize_lat, n_z), np.nan, dtype=float)
    half_c = tilesize_lon // 2 + 1
    half_a = tilesize_lat // 2 + 1
    tile[:half_c, :half_a, :] = block
    tile[half_c:, :half_a, :] = block[-2:0:-1, :, :]
    tile[:, half_a:tilesize_lat, :] = tile[:, tilesize_lat // 2 : 1 : -1, :]
    if seal_seam:
        tile = seal_periodic_cross_seam(tile)
    return tile


def extend_and_truncate(
    qf_tile: np.ndarray,
    n2mean: np.ndarray,
    depth: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extend forcing and $\\overline{N^2}$ below the GEM base, then average every three
    levels for a coarser vertical grid used in the omega inversion.
    """
    oo = N_DEPTH
    pp = oo + oo // 2
    n_cross, n_along, _ = qf_tile.shape
    qf_ext = np.full((n_cross, n_along, pp), np.nan, dtype=float)
    qf_ext[:, :, :oo] = qf_tile
    qf_ext[:, :, oo:] = qf_tile[:, :, oo - 1 : oo]

    n2_ext = np.full(pp, np.nan, dtype=float)
    n2_ext[:oo] = n2mean
    deep = (depth > 900) & (depth < 950)
    n2_ext[oo:] = np.nanmean(n2mean[deep])

    depth_ext = np.full(pp, np.nan, dtype=float)
    depth_ext[:oo] = depth
    depth_ext[oo:] = np.arange(991, 1480, 2)

    n_trunc = pp // 3
    qf_trunc = np.full((n_cross, n_along, n_trunc), np.nan, dtype=float)
    n2_trunc = np.full(n_trunc, np.nan, dtype=float)
    for i in range(n_trunc):
        sl = slice(i * 3, i * 3 + 3)
        block = qf_ext[:, :, sl]
        if np.isfinite(block).any():
            qf_trunc[:, :, i] = np.nanmean(block, axis=-1)
        n2_block = n2_ext[sl]
        if np.isfinite(n2_block).any():
            n2_trunc[i] = np.nanmean(n2_block)
    # MATLAB: depth_trunc = depth_extend(2:3:end-1)
    depth_trunc = depth_ext[1 : pp - 1 : 3]
    if depth_trunc.size > n_trunc:
        depth_trunc = depth_trunc[:n_trunc]
    elif depth_trunc.size < n_trunc:
        depth_trunc = depth_ext[1:pp:3][:n_trunc]
    return qf_trunc, n2_trunc, depth_trunc
