"""Default filesystem layout under the project root."""

from __future__ import annotations

from pathlib import Path

REPRO_ROOT = Path(__file__).resolve().parents[2]
VENDOR = REPRO_ROOT / "vendor" / "zenodo_19981347"
DATA = REPRO_ROOT / "data"
SWOT_DIR = DATA / "swot" / "_peek"
if not SWOT_DIR.exists():
    SWOT_DIR = DATA / "swot"
GEM_MAT = VENDOR / "sshbindens_681_final_curve_090.mat"
# Published FOCUS_data_trunc_*.nc buoyancy matches ``098`` (corr~1), not ``098_v2``.
# CreateDensityFields.m currently selects ``098_v2``; use the table that reproduces the deposit.
GEM_VAR = "sshbindens_681_final_curve_098"
OUT_DIR = DATA / "derived" / "python"
