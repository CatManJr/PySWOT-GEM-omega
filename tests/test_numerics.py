"""Small checks that the ported operators match the MATLAB definitions."""

from __future__ import annotations

import numpy as np

from swot_qg_omega.dynamics import project_uv_to_swath
from swot_qg_omega.io.gem import GemTable
from swot_qg_omega.io.swot import infer_orientation
from swot_qg_omega.qvector import q_vector
from swot_qg_omega.smooth import smooth_gaussian_2d
from swot_qg_omega.tile import mirror_tile_horizontal, seal_periodic_cross_seam


def test_geographic_projection_is_not_the_printed_matlab_pairing():
    ang = np.array([[np.nan, np.nan], [np.nan, 30.0]])
    u = np.array([[0.0, 0.0], [0.0, 1.0]])
    v = np.array([[0.0, 0.0], [0.0, 0.0]])
    cross, along = project_uv_to_swath(u, v, ang)
    rad = np.deg2rad(30.0)
    assert np.isclose(cross[1, 1], np.sin(rad))
    assert np.isclose(along[1, 1], -np.cos(rad))
    printed_cross = np.cos(rad)
    assert not np.isclose(cross[1, 1], printed_cross)


def test_q_vector_matches_create_density_fields():
    across_x = np.array([[[2.0]]])
    across_y = np.array([[[3.0]]])
    along_x = np.array([[[5.0]]])
    along_y = np.array([[[7.0]]])
    m2_across = np.array([[[11.0]]])
    m2_along = np.array([[[13.0]]])
    q1, q2 = q_vector(across_x, across_y, along_x, along_y, m2_across, m2_along)
    assert q1[0, 0, 0] == -(2.0 * 11.0 + 5.0 * 13.0)
    assert q2[0, 0, 0] == -(3.0 * 11.0 + 7.0 * 13.0)


def test_tile_seam_averages_first_and_last_cross_track_row():
    block = np.arange(4 * 3 * 2, dtype=float).reshape(4, 3, 2)
    tile = mirror_tile_horizontal(block, seal_seam=False)
    sealed = seal_periodic_cross_seam(tile)
    assert np.allclose(sealed[0], 0.5 * (tile[1] + tile[-1]))


def test_smoothdata_window_renormalizes_over_nans():
    field = np.array([[1.0, np.nan, 3.0], [1.0, 1.0, 1.0]], dtype=float)
    out = smooth_gaussian_2d(field, window=1)
    assert np.isclose(out[0, 0], 1.0)
    assert np.isnan(out[0, 1])
    assert np.isclose(out[0, 2], 3.0)


def test_gem_lookup_interpolates_between_bins():
    dens = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    table = GemTable(dens=dens, mid_ssh=np.array([-0.02, 0.0, 0.02]))
    mid = table.lookup_column(0.0)
    assert np.allclose(mid, dens[1])
    between = table.lookup_column(-0.01)
    assert np.allclose(between, 0.5 * dens[1] + 0.5 * dens[0])


def test_orientation_from_along_track_diagonal():
    lat = np.array([[-56.0, -55.0, -54.0]])
    lon = np.array([[150.0, 151.0, 152.0]])
    assert infer_orientation(lat, lon) == 1
    lat_se = np.array([[-54.0, -55.0, -56.0]])
    assert infer_orientation(lat_se, lon) == 2
