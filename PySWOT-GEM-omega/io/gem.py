"""Load a local GEM table mapping absolute dynamic topography to $\sigma_0(z)$."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from swot_qg_omega.constants import SSH_BIN_MAX, SSH_BIN_MIN, SSH_BIN_WIDTH
from swot_qg_omega.io.paths import GEM_MAT, GEM_VAR


@dataclass(frozen=True)
class GemTable:
    """GEM lookup dens[ssh_bin, depth]."""

    dens: np.ndarray
    mid_ssh: np.ndarray

    def lookup_column(self, adt: float) -> np.ndarray:
        """Linear interpolation between adjacent SSH bins."""
        dssh = self.mid_ssh[1] - self.mid_ssh[0]
        kk = int(np.searchsorted(self.mid_ssh, adt, side="right"))
        if kk <= 0:
            return self.dens[0].copy()
        if kk >= self.mid_ssh.size:
            return self.dens[-1].copy()
        sep = (self.mid_ssh[kk] - adt) / dssh
        return (1.0 - sep) * self.dens[kk] + sep * self.dens[kk - 1]


def ssh_bin_centers(
    ssh_min: float = SSH_BIN_MIN,
    ssh_max: float = SSH_BIN_MAX,
    width: float = SSH_BIN_WIDTH,
) -> np.ndarray:
    edges = np.arange(ssh_min, ssh_max + width * 0.5, width)
    return 0.5 * (edges[1:] + edges[:-1])


def load_gem(path: Path | None = None, var: str = GEM_VAR) -> GemTable:
    path = path or GEM_MAT
    raw = np.asarray(loadmat(path, squeeze_me=True)[var], dtype=float)
    # Stored arrays may be (depth, ssh) or (ssh, depth); normalize to (ssh, depth).
    dens = raw.T if raw.shape[0] > raw.shape[1] else raw
    mid_ssh = ssh_bin_centers()
    if dens.shape[0] != mid_ssh.size:
        raise ValueError(
            f"GEM ssh bins {dens.shape[0]} != expected mid_ssh {mid_ssh.size}"
        )
    return GemTable(dens=dens, mid_ssh=mid_ssh)
