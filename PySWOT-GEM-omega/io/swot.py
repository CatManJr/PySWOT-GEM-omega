"""Ingest regionally cropped SWOT L3 Expert KaRIn fields."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from swot_qg_omega.io.paths import SWOT_DIR

EPOCH_2000 = datetime(2000, 1, 1)


@dataclass
class SwotSnapshot:
    pass_id: str
    cycle_index: int
    yearday_2024: float
    lon: np.ndarray  # (cross-track, along-track)
    lat: np.ndarray
    adt: np.ndarray
    ucgos: np.ndarray
    vcgos: np.ndarray


def infer_orientation(lat: np.ndarray, lon: np.ndarray) -> int:
    """
    SWOT track sense used by CreateDensityFields.m.

    Returns 1 for a SW--NE along-track diagonal and 2 for NW--SE.
    The sense is the sign of $\\Delta\\mathrm{lat}\\times\\Delta\\mathrm{lon}$
    between the first and last finite samples on the center cross-track row.
    """
    i = lat.shape[0] // 2
    finite = np.isfinite(lat[i]) & np.isfinite(lon[i])
    idx = np.flatnonzero(finite)
    if idx.size < 2:
        return 1
    dlat = float(lat[i, idx[-1]] - lat[i, idx[0]])
    dlon = float(lon[i, idx[-1]] - lon[i, idx[0]])
    if dlon > 180.0:
        dlon -= 360.0
    elif dlon < -180.0:
        dlon += 360.0
    if dlat * dlon >= 0.0:
        return 1
    return 2


def list_focus_pass_ids(directory: Path | None = None) -> list[str]:
    """Pass ids present as FOCUS-cropped SWOT L3 Expert files."""
    directory = directory or SWOT_DIR
    ids = []
    for path in sorted(directory.glob("SWOT_L3_LR_SSH_Expert_*_v1.0.2_fcs.nc")):
        stem = path.name
        # SWOT_L3_LR_SSH_Expert_{pass}_v1.0.2_fcs.nc
        pass_id = stem.split("_")[5]
        ids.append(pass_id)
    return ids


def swot_filename(pass_id: str, directory: Path | None = None) -> Path:
    directory = directory or SWOT_DIR
    name = f"SWOT_L3_LR_SSH_Expert_{pass_id}_v1.0.2_fcs.nc"
    path = directory / name
    if not path.exists():
        alt = directory.parent / name
        if alt.exists():
            return alt
        raise FileNotFoundError(f"missing {path}; place L3 Expert files under {directory}")
    return path


def _cycle_datetime(time_row) -> datetime | None:
    row = np.asarray(time_row)
    if row.dtype.kind == "M":
        valid = row[~np.isnat(row)]
        if valid.size == 0:
            return None
        return valid[valid.size // 2].astype("datetime64[s]").tolist()
    valid = row[np.isfinite(row.astype(float))]
    if valid.size == 0:
        return None
    return EPOCH_2000 + timedelta(seconds=float(np.nanmedian(valid.astype(float))))


def list_january_cycles(pass_id: str, directory: Path | None = None) -> list[tuple[int, datetime]]:
    path = swot_filename(pass_id, directory)
    ds = xr.open_dataset(path)
    time = ds["time"].values
    out = []
    for ic in range(time.shape[0]):
        dt = _cycle_datetime(time[ic])
        if dt is not None and dt.year == 2024 and dt.month == 1:
            out.append((ic, dt))
    ds.close()
    return out


def load_snapshot(
    pass_id: str,
    cycle_index: int,
    directory: Path | None = None,
) -> SwotSnapshot:
    """
    Return one cycle on the swath grid.

    Arrays are ordered (num_pixels, num_lines) so axis 0 is cross-track and
    axis 1 is along-track.
    """
    path = swot_filename(pass_id, directory)
    ds = xr.open_dataset(path)
    lon = np.asarray(ds["longitude"].values.T, dtype=float)
    lat = np.asarray(ds["latitude"].values.T, dtype=float)
    ssha = np.asarray(ds["ssha"].values[cycle_index].T, dtype=float)
    mdt = np.asarray(ds["mdt"].values[cycle_index].T, dtype=float)
    ucgos = np.asarray(ds["ucgos_fit_13"].values[cycle_index].T, dtype=float)
    vcgos = np.asarray(ds["vcgos_fit_13"].values[cycle_index].T, dtype=float)
    dt = _cycle_datetime(ds["time"].values[cycle_index])
    ds.close()
    if dt is None:
        raise ValueError(f"no valid time for pass {pass_id} cycle {cycle_index}")
    yearday = (dt - datetime(2024, 1, 1)).total_seconds() / 86400.0

    return SwotSnapshot(
        pass_id=pass_id,
        cycle_index=cycle_index,
        yearday_2024=yearday,
        lon=lon,
        lat=lat,
        adt=ssha + mdt,
        ucgos=ucgos,
        vcgos=vcgos,
    )
