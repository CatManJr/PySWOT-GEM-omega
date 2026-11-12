"""Run local GEM + QG-omega reconstruction for one SWOT L3 snapshot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from swot_qg_omega.io.paths import OUT_DIR
from swot_qg_omega.io.swot import list_january_cycles
from swot_qg_omega.pipeline import run_pass_cycle, save_netcdf


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pass", dest="pass_id", default="155", help="SWOT pass id, e.g. 155")
    parser.add_argument(
        "--cycle",
        type=int,
        default=None,
        help="0-based cycle index; default = last January 2024 cycle for the pass",
    )
    parser.add_argument(
        "--orientation",
        type=int,
        choices=(1, 2),
        default=None,
        help="1 = SW-NE, 2 = NW-SE; default = infer from the along-track diagonal",
    )
    parser.add_argument("--no-w", action="store_true", help="Skip omega inversion (Q-vector only)")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    if args.cycle is None:
        jan = list_january_cycles(args.pass_id)
        if not jan:
            raise SystemExit(f"no January 2024 cycles for pass {args.pass_id}")
        args.cycle = jan[-1][0]
        print(f"Using January cycle index {args.cycle} ({jan[-1][1].isoformat()})")

    result = run_pass_cycle(
        args.pass_id,
        args.cycle,
        orientation=args.orientation,
        invert_w=not args.no_w,
    )
    out = args.out or (
        OUT_DIR / f"qg_omega_swath_{args.pass_id}_{result.cycle_index + 1}.nc"
    )
    save_netcdf(result, out)
    print(f"Wrote {out}")
    print(
        f"shapes buoyancy={result.buoyancy.shape} "
        f"qf_left={result.qf_left.shape} "
        f"w={'None' if result.w is None else result.w.shape}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
