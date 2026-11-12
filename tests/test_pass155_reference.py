"""Correctness against the one published FOCUS example (pass 155)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from swot_qg_omega.io.paths import OUT_DIR, VENDOR
from swot_qg_omega.io.swot import list_january_cycles
from swot_qg_omega.pipeline import run_pass_cycle
from swot_qg_omega.validate import overlap_corr, to_cross_along_depth

REFERENCE = VENDOR / "FOCUS_data_trunc_swath_155_9.nc"
SAVED = OUT_DIR / "qg_omega_swath_155_9.nc"
DIAGNOSTIC_MIN = {
    "buoyancy": 0.999,
    "vorticity": 0.999,
    "strain": 0.999,
    "divergence": 0.999,
    "N2mean_trunc": 0.999,
    "QF_trunc_left": 0.90,
    "QF_trunc_right": 0.90,
}


pytestmark = pytest.mark.skipif(
    not REFERENCE.exists(),
    reason="published FOCUS_data_trunc_swath_155_9.nc is not in vendor/",
)


def _score(reference: xr.Dataset, estimate: xr.Dataset, name: str) -> float:
    if name not in reference or name not in estimate:
        return float("nan")
    left = to_cross_along_depth(np.asarray(reference[name].values))
    right = np.asarray(estimate[name].values, dtype=float)
    if name in ("lon_buoy", "lat_buoy") and left.shape != right.shape and left.T.shape == right.shape:
        left = left.T
    if left.shape != right.shape:
        return float("nan")
    return overlap_corr(left, right)


def test_pass155_diagnostics_match_published_netcdf():
    cycles = list_january_cycles("155")
    if not cycles:
        pytest.skip("FOCUS SWOT L3 for pass 155 is not on disk")
    cycle_index = cycles[-1][0]
    result_path = OUT_DIR / "_pytest_pass155_now.nc"
    from swot_qg_omega.pipeline import save_netcdf

    result = run_pass_cycle("155", cycle_index, invert_w=False)
    save_netcdf(result, result_path)
    try:
        reference = xr.open_dataset(REFERENCE)
        estimate = xr.open_dataset(result_path)
        scores = {name: _score(reference, estimate, name) for name in DIAGNOSTIC_MIN}
        reference.close()
        estimate.close()
    finally:
        result_path.unlink(missing_ok=True)
    for name, floor in DIAGNOSTIC_MIN.items():
        assert scores[name] >= floor, f"{name} corr {scores[name]} < {floor}"


def test_saved_pass155_w_if_present():
    if not SAVED.exists():
        pytest.skip("run scripts/run_qg_omega.py --pass 155 to archive w")
    reference = xr.open_dataset(REFERENCE)
    estimate = xr.open_dataset(SAVED)
    if "w_vel" not in estimate:
        reference.close()
        estimate.close()
        pytest.skip("archived file has no w_vel")
    score = _score(reference, estimate, "w_vel")
    reference.close()
    estimate.close()
    assert score >= 0.99
