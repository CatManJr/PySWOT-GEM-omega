"""End-to-end local GEM, Q-vector forcing, and optional QG omega inversion."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from swot_qg_omega import density, dynamics, qvector, tile
from swot_qg_omega.constants import LAT_SUBSET_MAX, LAT_SUBSET_MIN
from swot_qg_omega.io.gem import GemTable, load_gem
from swot_qg_omega.io.swot import SwotSnapshot, infer_orientation, load_snapshot
from swot_qg_omega.omega import invert_tiled_omega, until_swath


@dataclass
class PipelineResult:
    pass_id: str
    cycle_index: int
    yearday_2024: float
    orientation: int
    lon: np.ndarray
    lat: np.ndarray
    adt: np.ndarray
    buoyancy: np.ndarray
    n2mean: np.ndarray
    vorticity: np.ndarray
    strain: np.ndarray
    divergence: np.ndarray
    qf: np.ndarray
    qf_left: np.ndarray
    qf_right: np.ndarray
    n2_trunc: np.ndarray
    depth_trunc: np.ndarray
    w: np.ndarray | None


def run_snapshot(
    snap: SwotSnapshot,
    gem: GemTable | None = None,
    *,
    orientation: int = 1,
    invert_w: bool = True,
) -> PipelineResult:
    gem = gem or load_gem()
    dens = density.reconstruct_sigma0(snap.adt, gem)
    lon, lat, dens, ucgos, vcgos = density.subset_latitude(
        snap.lon, snap.lat, dens, snap.ucgos, snap.vcgos
    )
    # subset ADT with the same along-track mask
    along_mask = (snap.lat[0, :] >= LAT_SUBSET_MIN) & (snap.lat[0, :] < LAT_SUBSET_MAX)
    adt = snap.adt[:, along_mask]
    lon, lat, dens, ucgos, vcgos = density.maybe_flip_orientation(
        lon, lat, dens, ucgos, vcgos, orientation
    )
    if orientation == 2:
        adt = adt[::-1, :][:, ::-1]


    dens = density.fill_surface_bottom_nans(dens)
    dens = density.enforce_stable(dens)

    buoyancy, n2mean, _rho0, _rhomean = density.buoyancy_and_n2(dens)
    m2_across, m2_along = density.horizontal_buoyancy_gradients(buoyancy)
    m2_across_sm, m2_along_sm = density.smooth_m2(m2_across, m2_along, buoyancy)

    ang = dynamics.swath_angle(lat, lon)
    u_cross, v_along = dynamics.project_uv_to_swath(ucgos, vcgos, ang)
    f = dynamics.coriolis(lat)
    full_across, full_along = dynamics.thermal_wind_velocity(
        m2_across_sm, m2_along_sm, f, u_cross, v_along
    )
    ax, ay, lx, ly = dynamics.velocity_gradients(full_across, full_along)
    ax, ay, lx, ly = dynamics.smooth_gradients(ax, ay, lx, ly)
    vorticity, strain, divergence = dynamics.vorticity_strain(ax, ay, lx, ly)

    q1, q2 = qvector.q_vector(ax, ay, lx, ly, m2_across_sm, m2_along_sm)
    qf = qvector.omega_forcing(q1, q2, buoyancy)

    edges = tile.swath_edges(orientation)
    i0, i1 = edges["left"]
    j0, j1 = edges["right"]
    left_block = qf[i0:i1, :, :]
    right_block = qf[j0:j1, :, :]
    qf_tile_l = tile.mirror_tile_horizontal(left_block)
    qf_tile_r = tile.mirror_tile_horizontal(right_block)

    depth = density.depth_coordinate()
    qf_left, n2_trunc, depth_trunc = tile.extend_and_truncate(qf_tile_l, n2mean, depth)
    qf_right, _, _ = tile.extend_and_truncate(qf_tile_r, n2mean, depth)

    w = None
    if invert_w:
        f0 = float(np.nanmean(f))
        w_left = invert_tiled_omega(qf_left, n2_trunc, depth_trunc, f0)
        w_right = invert_tiled_omega(qf_right, n2_trunc, depth_trunc, f0)
        w = until_swath(w_left, w_right, lon.shape[0], lon.shape[1], orientation)

    return PipelineResult(
        pass_id=snap.pass_id,
        cycle_index=snap.cycle_index,
        yearday_2024=snap.yearday_2024,
        orientation=orientation,
        lon=lon,
        lat=lat,
        adt=adt,
        buoyancy=buoyancy,
        n2mean=n2mean,
        vorticity=vorticity,
        strain=strain,
        divergence=divergence,
        qf=qf,
        qf_left=qf_left,
        qf_right=qf_right,
        n2_trunc=n2_trunc,
        depth_trunc=depth_trunc,
        w=w,
    )


def run_pass_cycle(
    pass_id: str,
    cycle_index: int,
    *,
    orientation: int | None = None,
    invert_w: bool = True,
) -> PipelineResult:
    snap = load_snapshot(pass_id, cycle_index)
    if orientation is None:
        orientation = infer_orientation(snap.lat, snap.lon)
    return run_snapshot(snap, orientation=orientation, invert_w=invert_w)


def to_dataset(result: PipelineResult) -> xr.Dataset:
    data_vars = {
        "buoyancy": (("cross", "along", "depth"), result.buoyancy),
        "vorticity": (("cross", "along", "depth"), result.vorticity),
        "strain": (("cross", "along", "depth"), result.strain),
        "divergence": (("cross", "along", "depth"), result.divergence),
        "qf": (("cross", "along", "depth"), result.qf),
        "QF_trunc_left": (("cross_tile_l", "along_tile", "depth_trunc"), result.qf_left),
        "QF_trunc_right": (("cross_tile_r", "along_tile", "depth_trunc"), result.qf_right),
        "N2mean_trunc": (("depth_trunc",), result.n2_trunc),
        "lon_buoy": (("cross", "along"), result.lon),
        "lat_buoy": (("cross", "along"), result.lat),
        "depth_trunc": (("depth_trunc",), result.depth_trunc),
        "YearDay": (("yearday",), np.array([result.yearday_2024])),
        "orientation": (("swathindex",), np.array([result.orientation], dtype=float)),
    }
    if result.w is not None:
        data_vars["w_vel"] = (("cross", "along", "depth_trunc"), result.w)
    return xr.Dataset(
        data_vars,
        coords={
            "depth": density.depth_coordinate(),
        },
        attrs={
            "pass_id": result.pass_id,
            "cycle_index": result.cycle_index,
            "pipeline": "swot_qg_omega: Python port of CreateDensityFields.m + xinvert omega",
        },
    )


def save_netcdf(result: PipelineResult, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    to_dataset(result).to_netcdf(path)
