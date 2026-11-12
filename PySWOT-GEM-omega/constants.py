"""Physical constants and default numerical parameters for the QG-omega pipeline."""

from __future__ import annotations

import numpy as np

G = 9.81
OMEGA_EARTH = 7.292115e-5
GRID_SPACING_M = 2000.0

# Local GEM (SSH -> sigma0) binning
SSH_BIN_MIN = -1.10
SSH_BIN_MAX = 0.20
SSH_BIN_WIDTH = 0.02

# Vertical grid for reconstructed density (m)
DEPTH_M = np.arange(11, 990, 2, dtype=float)
N_DEPTH = DEPTH_M.size

# MATLAB smoothdata(..., 'gaussian', 8): sigma = window/5, truncated window
GAUSS_WINDOW = 8
GAUSS_SIGMA = GAUSS_WINDOW / 5.0  # documented MATLAB ratio; used only as reference

# Along-track latitude band retained for the ACC FOCUS box
LAT_SUBSET_MIN = -58.0
LAT_SUBSET_MAX = -52.0

# Cross-track indices of left/right KaRIn swaths (1-based inclusive; convert at use)
SWATH_EDGES_ORIENT1 = {
    "left": (5, 31),
    "right": (40, 66),
}
SWATH_EDGES_ORIENT2 = {
    "left": (5, 30),
    "right": (40, 65),
}
