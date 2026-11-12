"""Compare Python QG-omega output against the published MATLAB reference NetCDF."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

REF_DEFAULT = ROOT / "vendor" / "zenodo_19981347" / "FOCUS_data_trunc_swath_155_9.nc"
PY_DEFAULT = ROOT / "data" / "derived" / "python" / "qg_omega_swath_155_9.nc"


def _to_cross_along_depth(arr: np.ndarray, name: str) -> np.ndarray:
    """Normalize MATLAB NetCDF (depth, along, cross) to (cross, along, depth)."""
    a = np.asarray(arr, dtype=float)
    if a.ndim != 3:
        return a
    # Published file uses (depth, along, cross) for full-swath fields.
    if a.shape[0] in (245, 490) and a.shape[-1] in (52, 69):
        return np.transpose(a, (2, 1, 0))
    return a


def _stats(ref: np.ndarray, py: np.ndarray, label: str) -> dict:
    mask = np.isfinite(ref) & np.isfinite(py)
    n = int(mask.sum())
    if n < 10:
        print(f"{label:18s}  insufficient overlap (n={n})")
        return {"label": label, "n": n, "corr": np.nan, "rmse": np.nan, "bias": np.nan}
    r = ref[mask]
    p = py[mask]
    corr = float(np.corrcoef(r, p)[0, 1]) if r.std() > 0 and p.std() > 0 else np.nan
    rmse = float(np.sqrt(np.mean((p - r) ** 2)))
    bias = float(np.mean(p - r))
    rel = rmse / (float(np.std(r)) + 1e-30)
    print(
        f"{label:18s}  n={n:9d}  corr={corr:7.4f}  rmse={rmse:.4e}  "
        f"bias={bias:.4e}  rmse/std_ref={rel:.3f}  "
        f"|ref|_mean={np.mean(np.abs(r)):.4e}  |py|_mean={np.mean(np.abs(p)):.4e}"
    )
    return {
        "label": label,
        "n": n,
        "corr": corr,
        "rmse": rmse,
        "bias": bias,
        "rmse_over_std": rel,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", type=Path, default=REF_DEFAULT)
    parser.add_argument("--python", dest="py_path", type=Path, default=None)
    args = parser.parse_args()

    candidates = [
        args.py_path,
        ROOT / "data" / "derived" / "python" / "qg_omega_swath_155_9.nc",
        ROOT / "data" / "derived" / "python" / "FOCUS_python_swath_155_9.nc",
    ]
    py_path = next((p for p in candidates if p is not None and p.exists()), None)
    if py_path is None:
        raise SystemExit("no Python output NetCDF found; run scripts/run_qg_omega.py first")
    if not args.ref.exists():
        raise SystemExit(f"missing MATLAB reference: {args.ref}")

    print(f"MATLAB ref: {args.ref}")
    print(f"Python out: {py_path}")
    ref = xr.open_dataset(args.ref)
    py = xr.open_dataset(py_path)

    print("\n=== field comparison (Python vs MATLAB) ===")
    pairs = [
        ("buoyancy", "buoyancy"),
        ("vorticity", "vorticity"),
        ("strain", "strain"),
        ("divergence", "divergence"),
        ("QF_trunc_left", "QF_trunc_left"),
        ("QF_trunc_right", "QF_trunc_right"),
        ("N2mean_trunc", "N2mean_trunc"),
        ("w_vel", "w_vel"),
        ("lon_buoy", "lon_buoy"),
        ("lat_buoy", "lat_buoy"),
    ]
    results = []
    for rname, pname in pairs:
        if rname not in ref or pname not in py:
            print(f"{rname:18s}  skipped (missing in one file)")
            continue
        r = _to_cross_along_depth(ref[rname].values, rname)
        p = np.asarray(py[pname].values, dtype=float)
        # lon/lat in ref are (along, cross) after netcdf read of (cross, along) write
        if rname in ("lon_buoy", "lat_buoy") and r.ndim == 2:
            if r.shape != p.shape and r.T.shape == p.shape:
                r = r.T
        if r.shape != p.shape:
            print(f"{rname:18s}  shape mismatch ref={r.shape} py={p.shape}")
            continue
        results.append(_stats(r, p, rname))

    print("\n=== scalar metadata ===")
    for key in ("YearDay", "orientation"):
        if key in ref and key in py:
            print(f"{key}: ref={np.asarray(ref[key].values).ravel()}  py={np.asarray(py[key].values).ravel()}")

    # Summary gate for the fields that should match well
    key = {r["label"]: r for r in results}
    buoy = key.get("buoyancy", {})
    print("\n=== verdict ===")
    if buoy.get("corr", 0) >= 0.95:
        print("Density/buoyancy stage: CLOSE (corr >= 0.95)")
    else:
        print("Density/buoyancy stage: DIVERGES from MATLAB reference")
    if "w_vel" in key and np.isfinite(key["w_vel"].get("corr", np.nan)):
        print(f"Omega w stage: corr={key['w_vel']['corr']:.4f}")
    elif "w_vel" not in key:
        print("Omega w stage: not present in Python file")

    ref.close()
    py.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
