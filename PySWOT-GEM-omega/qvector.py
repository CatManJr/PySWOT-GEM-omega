"""Q-vector construction and omega-equation forcing (CreateDensityFields.m)."""

from __future__ import annotations

import numpy as np

from swot_qg_omega.constants import GRID_SPACING_M
from swot_qg_omega.smooth import smooth_gaussian_2d


def q_vector(
    across_x: np.ndarray,
    across_y: np.ndarray,
    along_x: np.ndarray,
    along_y: np.ndarray,
    m2_across: np.ndarray,
    m2_along: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Hoskins Q-vector on the swath grid:

    ``Q1 = -u_x M2_x - v_x M2_y``, ``Q2 = -u_y M2_x - v_y M2_y``
    with $(u,v)=(\\mathrm{across},\\mathrm{along})$.
    """
    q1 = -across_x * m2_across - along_x * m2_along
    q2 = -across_y * m2_across - along_y * m2_along
    return q1, q2


def divergence_of_q(q1: np.ndarray, q2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centered $\\partial_x Q_1$ and $\\partial_y Q_2$ on the 2 km grid."""
    dx = GRID_SPACING_M
    q1x = np.full_like(q1, np.nan)
    q2y = np.full_like(q2, np.nan)
    q1x[1:-1, :, :] = (q1[2:, :, :] - q1[:-2, :, :]) / (2.0 * dx)
    q2y[:, 1:-1, :] = (q2[:, 2:, :] - q2[:, :-2, :]) / (2.0 * dx)
    return q1x, q2y


def omega_forcing(q1: np.ndarray, q2: np.ndarray, buoyancy: np.ndarray) -> np.ndarray:
    """
    ``QFsmooth = 2*(Q1_x_sm + Q2_y_sm)`` as in CreateDensityFields.m.

    Nadir-gap cross-track indices are restored from
    ``isnan(buoyancy(:,200,200))`` (MATLAB 1-based).
    """
    q1x, q2y = divergence_of_q(q1, q2)
    q1x_sm = np.empty_like(q1x)
    q2y_sm = np.empty_like(q2y)
    for k in range(q1x.shape[-1]):
        q1x_sm[..., k] = smooth_gaussian_2d(q1x[..., k])
        q2y_sm[..., k] = smooth_gaussian_2d(q2y[..., k])
    qf = 2.0 * (q1x_sm + q2y_sm)
    # MATLAB: kk = find(isnan(buoyancy(:,200,200))); QFsmooth(kk,:,:) = NaN
    along_i = min(199, buoyancy.shape[1] - 1)
    depth_i = min(199, buoyancy.shape[2] - 1)
    gap = np.isnan(buoyancy[:, along_i, depth_i])
    qf[gap, :, :] = np.nan
    return qf
