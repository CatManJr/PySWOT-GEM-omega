"""Local GEM density reconstruction, buoyancy, and stratification."""

from __future__ import annotations

import numpy as np

from swot_qg_omega.constants import (
    DEPTH_M,
    G,
    GRID_SPACING_M,
    LAT_SUBSET_MAX,
    LAT_SUBSET_MIN,
    N_DEPTH,
)
from swot_qg_omega.io.gem import GemTable
from swot_qg_omega.smooth import smooth_gaussian_2d


def reconstruct_sigma0(adt: np.ndarray, gem: GemTable) -> np.ndarray:
    """Map ADT to $\sigma_0(z)$. Returns (cross, along, depth)."""
    ny, nx = adt.shape
    out = np.full((ny, nx, N_DEPTH), np.nan, dtype=float)
    for i in range(ny):
        for j in range(nx):
            if np.isfinite(adt[i, j]):
                out[i, j, :] = gem.lookup_column(float(adt[i, j]))
    return out


def subset_latitude(
    lon: np.ndarray,
    lat: np.ndarray,
    dens: np.ndarray,
    ucgos: np.ndarray,
    vcgos: np.ndarray,
    lat_min: float = LAT_SUBSET_MIN,
    lat_max: float = LAT_SUBSET_MAX,
) -> tuple[np.ndarray, ...]:
    """Retain the along-track band $\mathrm{lat\_min} \le \mathrm{lat} < \mathrm{lat\_max}$."""
    along_mask = (lat[0, :] >= lat_min) & (lat[0, :] < lat_max)
    return (
        lon[:, along_mask],
        lat[:, along_mask],
        dens[:, along_mask, :],
        ucgos[:, along_mask],
        vcgos[:, along_mask],
    )


def maybe_flip_orientation(
    lon: np.ndarray,
    lat: np.ndarray,
    dens: np.ndarray,
    ucgos: np.ndarray,
    vcgos: np.ndarray,
    orientation: int,
) -> tuple[np.ndarray, ...]:
    """
    Orientation flag: 1 keeps SW--NE ordering; 2 reverses both axes for NW--SE tracks
    so left/right KaRIn edges share a common index convention.
    """
    if orientation != 2:
        return lon, lat, dens, ucgos, vcgos
    lon = lon[::-1, :][:, ::-1]
    lat = lat[::-1, :][:, ::-1]
    dens = dens[::-1, :, :][:, ::-1, :]
    ucgos = ucgos[::-1, :][:, ::-1]
    vcgos = vcgos[::-1, :][:, ::-1]
    return lon, lat, dens, ucgos, vcgos


def fill_surface_bottom_nans(dens: np.ndarray) -> np.ndarray:
    """Fill near-surface gaps by nearest deeper value; extrapolate near the bottom."""
    out = dens.copy()
    ny, nx, nz = out.shape
    for i in range(ny):
        for j in range(nx):
            col = out[i, j, :]
            top = np.where(np.isnan(col[:30]))[0]
            if top.size:
                maxkk = int(top.max())
                if maxkk + 1 < nz and np.isfinite(col[maxkk + 1]):
                    col[: maxkk + 1] = col[maxkk + 1]
            bot = np.where(np.isnan(col[400:490]))[0]
            if bot.size:
                minkk = int(bot.min()) + 400
                if minkk > 401 and np.isfinite(col[minkk - 1]):
                    x = np.arange(401, minkk)
                    y = col[401:minkk]
                    if np.sum(np.isfinite(y)) >= 2:
                        coef = np.polyfit(x[np.isfinite(y)], y[np.isfinite(y)], 1)
                        col[minkk:490] = np.polyval(coef, np.arange(minkk, 490))
                    else:
                        col[minkk:490] = np.nan
                else:
                    col[minkk:490] = np.nan
            out[i, j, :] = col
    return out


def enforce_stable(dens: np.ndarray) -> np.ndarray:
    """Sort each column so potential density increases with depth."""
    return np.sort(dens, axis=-1)


def buoyancy_and_n2(dens: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Convert $\sigma_0$ to buoyancy $b$ and domain-mean $N^2(z)$.

    Uses $b = g(\rho_0 - \rho)/\rho_0$ with $\rho = 1000 + \sigma_0$.
    """
    rhomean = np.nanmean(dens, axis=(0, 1))
    rho0 = 1000.0 + float(np.nanmean(rhomean))
    buoyancy = (rho0 - (1000.0 + dens)) / rho0 * G
    n2mean = np.full(N_DEPTH, np.nan, dtype=float)
    n2mean[1:-1] = (rhomean[2:] - rhomean[:-2]) / 4.0 / rho0 * G
    return buoyancy, n2mean, rho0, rhomean


def horizontal_buoyancy_gradients(buoyancy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centered $M^2$ on the 2 km grid: across-track and along-track components."""
    dx = GRID_SPACING_M
    m2_across = np.full_like(buoyancy, np.nan)
    m2_along = np.full_like(buoyancy, np.nan)
    m2_across[1:-1, :, :] = (buoyancy[2:, :, :] - buoyancy[:-2, :, :]) / (2.0 * dx)
    m2_along[:, 1:-1, :] = (buoyancy[:, 2:, :] - buoyancy[:, :-2, :]) / (2.0 * dx)
    return m2_across, m2_along


def smooth_m2(
    m2_across: np.ndarray,
    m2_along: np.ndarray,
    buoyancy: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian-smooth $M^2$ and restore the nadir gap from buoyancy missing values."""
    m2a = np.empty_like(m2_across)
    m2l = np.empty_like(m2_along)
    for k in range(m2_across.shape[-1]):
        m2a[..., k] = smooth_gaussian_2d(m2_across[..., k])
        m2l[..., k] = smooth_gaussian_2d(m2_along[..., k])
    # MATLAB: kk = find(isnan(buoyancy(:,1,200)))
    gap = np.isnan(buoyancy[:, 0, 199])
    m2a[gap, :, :] = np.nan
    m2l[gap, :, :] = np.nan
    return m2a, m2l


def depth_coordinate() -> np.ndarray:
    return DEPTH_M.copy()
