# SWOT local GEM and QG omega

Python replacement for the MATLAB half of Thompson et al. (*Communications Earth & Environment*, 2026, doi:10.1038/s43247-026-03673-y). The quasi-geostrophic omega inversion for vertical velocity $w$ stays in [`xinvert`](https://github.com/miniufo/xinvert), which is the solver named in that paper.

The package reconstructs upper-ocean density on a SWOT KaRIn swath from a local gravest empirical mode (GEM), builds the thermal-wind velocity and the Hoskins Q-vector, tiles each swath for an elliptic solve, and writes $w$ back onto the swath. It is a migration of the deposited diagnostic script, not a new dynamical model and not a substitute for the published SWOT or glider observations.

## What moved, and what did not

Thompson et al. split the calculation in the same way this repository does. `CreateDensityFields.m` (Zenodo doi:10.5281/zenodo.19981347) reads SWOT L3 absolute dynamic topography and cyclogeostrophic surface velocity, looks up $\sigma_0(z)$ from a precomputed GEM, and forms buoyancy, $M^2$, thermal-wind shear, vorticity, strain, and the omega forcing $2\nabla\cdot\mathbf{Q}$. That script does not invert for $w$. It reads a separate file, written by `xinvert`, and stitches the left and right solutions onto the swath. Code availability in the paper states the same division.

This package therefore replaces only the MATLAB diagnostic path. `numpy`, `scipy`, `numba`, and `xarray` implement the GEM lookup, centered differences, Gaussian smoothing, Q-vector, and mirror tiling. `xinvert.invert_omega` remains the elliptic solver, with Dirichlet conditions $w=0$ at the top and bottom of the extended depth grid. Horizontal periodicity is requested on both tiled axes. Current `xinvert` kernels enforce a periodic condition only on the last axis; that limit is reported at runtime and is not hidden by a second solver.

Two lines in the deposited `.m` file do not reproduce the published example NetCDF, so the port follows the product rather than those lines.

The printed rotation `cross = u\cos\theta + v\sin\theta` anticorrelates with the published vorticity. Here $\theta$ is the `sw_dist` bearing of the cross-track chord from north, and $(u,v)$ is eastward and northward velocity. The across-track unit vector is $(\sin\theta,\cos\theta)$, so

$$u_\perp = u\sin\theta + v\cos\theta, \qquad v_\parallel = -u\cos\theta + v\sin\theta.$$

With that projection, pass 155 vorticity, strain, and divergence match the published fields.

The script selects GEM table `sshbindens_681_final_curve_098_v2`. Published buoyancy on the same pass matches `sshbindens_681_final_curve_098` in the same `.mat` file. This package loads `098`. Smoothing follows `smoothdata(..., 'gaussian', 8)`: a truncated Gaussian with $\sigma = \mathrm{window}/5$, an even window centered on the current and previous sample, NaNs omitted, and weights renormalized on the samples that remain.

On the one archived example (`FOCUS_data_trunc_swath_155_9.nc`, SWOT pass 155, late January 2024), correlations against this pipeline are $1.000$ for buoyancy, vorticity, strain, divergence, and $\overline{N^2}$; about $0.95$--$0.96$ for the tiled Q-vector forcing; and $0.995$ for $w$. The remaining forcing difference, and a low bias in $|w|$ of about one quarter even when the published forcing is passed through the same `xinvert` call, are consistent with the unpublished `xinvert` options used for the archive and with the missing second periodic boundary. They are not removed by rescaling.

The GEM is the campaign table. Glider profiles are not refit here. Only one vertical-velocity file is in the Zenodo deposit; other passes can be regenerated with the batch script below. SWOT swaths in this tree are the FOCUS crops shipped with that deposit, not a new AVISO download.

## Data products

| Role | Product | Reference |
| --- | --- | --- |
| Reference code, GEM, SWOT crops, example $w$ | Zenodo | https://doi.org/10.5281/zenodo.19981347 |
| Seagliders SG659 and SG681 | Zenodo | https://doi.org/10.5281/zenodo.17252347 |
| SWOT L3 KaRIn LR Expert v1.0.2 | AVISO/DUACS | https://doi.org/10.24400/527896/A01-2023.018 |
| DUACS L4 SLA/ADT | CMEMS/DUACS | https://doi.org/10.48670/moi-00148 |
| QG omega solver | `xinvert` | https://github.com/miniufo/xinvert |

## Layout

```
swot_qg_omega/
  swot_qg_omega/            # package
  tests/
  vendor/zenodo_19981347/   # deposited scripts and reference files
  data/swot/                # FOCUS SWOT L3 Expert swaths
  data/glider/              # Seaglider archives
  data/duacs/
  data/derived/python/      # pipeline outputs
  scripts/
  env/
```

## Setup

```bash
conda env create -f env/environment.yml
conda activate swot-qg-omega

# or
python -m pip install -r env/requirements.txt

python scripts/download_zenodo.py
python scripts/download_zenodo.py --full
```

## Run

```bash
python scripts/run_qg_omega.py --pass 155 --no-w
python scripts/run_qg_omega.py --pass 155
python scripts/compare_matlab_python.py
python scripts/run_focus_passes.py
python scripts/run_focus_passes.py --with-w
pytest
```

`--no-w` writes density, velocity diagnostics, and Q-vector forcing. The default single-pass command also calls `xinvert`, which takes on the order of 15 minutes per swath. Orientation defaults to the along-track diagonal (SW--NE or NW--SE). `run_focus_passes.py` uses the last January 2024 cycle of each FOCUS file under `data/swot/` and, unless `--with-w` is set, skips the inversion. `pytest` checks the operators and, when the reference NetCDF and pass 155 SWOT file are present, recomputes the diagnostic fields against that NetCDF. The archived $w$ correlation is checked only if `data/derived/python/qg_omega_swath_155_9.nc` already contains `w_vel`.
