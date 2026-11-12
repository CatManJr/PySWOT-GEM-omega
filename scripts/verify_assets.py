"""Verify local FOCUS SWOT / glider / GEM assets against Zenodo size and md5."""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "vendor" / "zenodo_19981347"
GLIDER_DIR = ROOT / "data" / "glider"
ALT_DIRS = [
    ROOT / "data" / "derived",
    ROOT / "data" / "duacs",
    ROOT / "data" / "swot",
]


def fetch_record(record_id: int) -> dict:
    url = f"https://zenodo.org/api/records/{record_id}"
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def md5_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def resolve_local(key: str, primary: Path) -> Path | None:
    candidate = primary / key
    if candidate.exists():
        return candidate
    for folder in ALT_DIRS:
        alt = folder / key
        if alt.exists():
            return alt
    return None


def check_record(record: dict, primary: Path, label: str, do_md5: bool) -> tuple[int, int, int]:
    print(f"=== {label} ({record.get('doi')}) ===")
    ok = fail = miss = 0
    for item in sorted(record["files"], key=lambda x: x["key"]):
        key = item["key"]
        expected_size = int(item["size"])
        checksum = item.get("checksum", "")
        expected_md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else None
        path = resolve_local(key, primary)
        if path is None:
            print(f"  MISS  {key}")
            miss += 1
            continue
        got_size = path.stat().st_size
        size_ok = got_size == expected_size
        md5_ok = True
        note = ""
        if do_md5 and expected_md5:
            print(f"  ....  hashing {key} ({got_size / 1e6:.1f} MB) ...", flush=True)
            got_md5 = md5_file(path)
            md5_ok = got_md5 == expected_md5
            if not md5_ok:
                note += f" md5={got_md5}/{expected_md5}"
        if not size_ok:
            note += f" size={got_size}/{expected_size}"
        status = "OK" if size_ok and md5_ok else "FAIL"
        if status == "OK":
            ok += 1
        else:
            fail += 1
        rel = path.relative_to(ROOT).as_posix()
        print(f"  {status:4}  {key}  {got_size / 1e6:.1f} MB  [{rel}]{note}")
    print(f"  summary: ok={ok} fail={fail} miss={miss}\n")
    return ok, fail, miss


def main() -> int:
    do_md5 = "--no-md5" not in sys.argv
    code = fetch_record(19981347)
    glider = fetch_record(17252347)
    a = check_record(code, CODE_DIR, "zenodo code/data", do_md5)
    b = check_record(glider, GLIDER_DIR, "zenodo glider", do_md5)

    print("=== layout mirrors ===")
    mirrors = {
        "duacs L4": ROOT / "data" / "duacs" / "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc",
        "mean SSH": ROOT / "data" / "duacs" / "means_SSH_1993_2017.nc",
        "GEM mat": ROOT / "data" / "derived" / "sshbindens_681_final_curve_090.mat",
        "example w": ROOT / "data" / "derived" / "FOCUS_data_trunc_swath_155_9.nc",
        "SWOT zip": ROOT / "data" / "swot" / "SWOT_Swaths_ThompsonetalCEE26.zip",
    }
    for name, path in mirrors.items():
        print(f"  {'OK' if path.exists() else 'MISS':4}  {name}: {path.relative_to(ROOT) if path.exists() else path}")

    print("\n=== python libs ===")
    for name in ["xinvert", "xarray", "gsw", "glidertools", "numba", "dask", "netCDF4"]:
        try:
            mod = __import__(name)
            print(f"  OK    {name} {getattr(mod, '__version__', '?')}")
        except Exception as exc:  # noqa: BLE001
            print(f"  FAIL  {name}: {exc}")

    bad = a[1] + a[2] + b[1] + b[2]
    print(f"\nTOTAL fail+miss: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
