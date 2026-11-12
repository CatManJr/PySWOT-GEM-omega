"""Cyclogeostrophic surface projection, thermal wind, and strain diagnostics."""

from __future__ import annotations

import numpy as np

from swot_qg_omega.constants import GRID_SPACING_M, OMEGA_EARTH
from swot_qg_omega.smooth import smooth_gaussian_2d


def coriolis(lat: np.ndarray) -> np.ndarray:
    return 2.0 * OMEGA_EARTH * np.sin(np.deg2rad(lat))


def swath_angle(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """
    Cross-track pair angle from seawater ``sw_dist`` (CreateDensityFields.m).

    ``[dist, ang] = sw_dist([lat(i-1,j) lat(i+1,j)], [lon(i-1,j) lon(i+1,j)], 'km')``
    with ``ang = atan2(dlon*cos(lat2), dlat)*180/pi``.
    """
    ang = np.full(lat.shape, np.nan, dtype=float)
    for i in range(1, lat.shape[0] - 1):
        for j in range(1, lat.shape[1] - 1):
            lat1, lon1 = lat[i - 1, j], lon[i - 1, j]
            lat2, lon2 = lat[i + 1, j], lon[i + 1, j]
            if not np.isfinite([lat1, lon1, lat2, lon2]).all():
                continue
            dlat = lat2 - lat1
            dlon = lon2 - lon1
            if dlon > 180.0:
                dlon -= 360.0
            elif dlon < -180.0:
                dlon += 360.0
            ang[i, j] = np.rad2deg(np.arctan2(dlon * np.cos(np.deg2rad(lat2)), dlat))
    return ang


def project_uv_to_swath(
    u: np.ndarray,
    v: np.ndarray,
    ang_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Project eastward/northward velocity onto across- and along-swath axes.

    ``sw_dist`` returns bearing $\theta$ of the cross-track chord from north, so
    the across-track unit vector in $(u,v)=(\mathrm{east},\mathrm{north})$ is
    $(\sin\theta,\cos\theta)$ and the along-track unit vector (90$^\circ$ CCW) is
    $(-\cos\theta,\sin\theta)$:

    ``cross = u*sin(ang) + v*cos(ang)``
    ``along = -u*cos(ang) + v*sin(ang)``

    This matches Thompson et al. (2026) across-/along-swath diagnostics and the
    published Zenodo vorticity. The trig pairing printed in CreateDensityFields.m
    (``u*cos+v*sin``, ``-u*sin+v*cos``) anticorrelates with that product.
    """
    rad = np.deg2rad(ang_deg)
    cross = u * np.sin(rad) + v * np.cos(rad)
    along = -u * np.cos(rad) + v * np.sin(rad)
    return cross, along


def thermal_wind_velocity(
    m2_across_sm: np.ndarray,
    m2_along_sm: np.ndarray,
    f: np.ndarray,
    u_cross_s: np.ndarray,
    v_along_s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Integrate thermal-wind shear from smoothed $M^2$ and add surface velocity.

    Vertical spacing is $2\,\mathrm{m}$ on the GEM depth grid.
    """
    ff = np.broadcast_to(f[:, :, None], m2_across_sm.shape)
    shear_along = -np.cumsum(m2_across_sm, axis=-1) * 2.0 / ff
    shear_across = np.cumsum(m2_along_sm, axis=-1) * 2.0 / ff
    full_across = shear_across + u_cross_s[:, :, None]
    full_along = shear_along + v_along_s[:, :, None]
    return full_across, full_along


def velocity_gradients(
    full_across: np.ndarray,
    full_along: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Centered horizontal derivatives of cross- and along-track velocity."""
    dx = GRID_SPACING_M
    ax = np.full_like(full_across, np.nan)
    ay = np.full_like(full_across, np.nan)
    lx = np.full_like(full_along, np.nan)
    ly = np.full_like(full_along, np.nan)
    ax[1:-1, :, :] = (full_across[2:, :, :] - full_across[:-2, :, :]) / (2.0 * dx)
    ay[:, 1:-1, :] = (full_across[:, 2:, :] - full_across[:, :-2, :]) / (2.0 * dx)
    lx[1:-1, :, :] = (full_along[2:, :, :] - full_along[:-2, :, :]) / (2.0 * dx)
    ly[:, 1:-1, :] = (full_along[:, 2:, :] - full_along[:, :-2, :]) / (2.0 * dx)
    return ax, ay, lx, ly


def smooth_gradients(
    ax: np.ndarray,
    ay: np.ndarray,
    lx: np.ndarray,
    ly: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    out = []
    for field in (ax, ay, lx, ly):
        sm = np.empty_like(field)
        for k in range(field.shape[-1]):
            sm[..., k] = smooth_gaussian_2d(field[..., k])
        out.append(sm)
    return tuple(out)  # type: ignore[return-value]


def vorticity_strain(
    across_x: np.ndarray,
    across_y: np.ndarray,
    along_x: np.ndarray,
    along_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vorticity = along_x - across_y
    strain_n = across_x - along_y
    strain_s = along_x + across_y
    strain = np.sqrt(strain_n**2 + strain_s**2)
    divergence = across_x + along_y
    return vorticity, strain, divergence
