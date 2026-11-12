"""Run the local GEM + QG-omega pipeline on every FOCUS SWOT pass on disk."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from swot_qg_omega.io.paths import OUT_DIR
from swot_qg_omega.io.swot import infer_orientation, list_focus_pass_ids, list_january_cycles, load_snapshot
from swot_qg_omega.pipeline import run_pass_cycle, save_netcdf


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--with-w",
        action="store_true",
        help="Invert w with xinvert (slow). Default writes diagnostics and Q only.",
    )
    parser.add_argument(
        "--passes",
        nargs="*",
        default=None,
        help="Pass ids to run; default = every FOCUS L3 file under data/swot",
    )
    args = parser.parse_args()

    pass_ids = args.passes or list_focus_pass_ids()
    if not pass_ids:
        raise SystemExit("no FOCUS SWOT L3 files found")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = OUT_DIR / "focus_batch_summary.txt"
    lines = ["pass_id cycle yearday orientation out finite_buoyancy"]
    for pass_id in pass_ids:
        cycles = list_january_cycles(pass_id)
        if not cycles:
            print(f"pass {pass_id}: no January 2024 cycle")
            continue
        cycle_index, when = cycles[-1]
        snap = load_snapshot(pass_id, cycle_index)
        orientation = infer_orientation(snap.lat, snap.lon)
        print(
            f"pass {pass_id} cycle {cycle_index} {when.isoformat()} "
            f"orientation {orientation} invert_w={args.with_w}"
        )
        result = run_pass_cycle(
            pass_id,
            cycle_index,
            orientation=orientation,
            invert_w=args.with_w,
        )
        stem = f"qg_omega_swath_{pass_id}_{result.cycle_index + 1}"
        # Diagnostics-only runs must not replace an archived xinvert solution.
        out = OUT_DIR / (f"{stem}.nc" if args.with_w else f"{stem}_diag.nc")
        save_netcdf(result, out)
        finite = float(np.isfinite(result.buoyancy).mean())
        line = (
            f"{pass_id} {result.cycle_index + 1} {result.yearday_2024:.4f} "
            f"{result.orientation} {out.name} {finite:.4f}"
        )
        lines.append(line)
        print(f"  wrote {out}")

    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Summary {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
