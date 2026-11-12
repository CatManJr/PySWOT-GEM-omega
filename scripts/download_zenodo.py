"""Download published FOCUS SWOT / glider / GEM assets from Zenodo.

Records (paper Data/Code availability):
  - Reference code, GEM table, FOCUS SWOT package: https://doi.org/10.5281/zenodo.19981347
  - Seagliders SG659 / SG681: https://doi.org/10.5281/zenodo.17252347

Remote filenames below are unchanged archive keys.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "DATA_MANIFEST.json"

ZENODO_CODE = 19981347
ZENODO_GLIDER = 17252347

# Prefer method inputs first; large archives optional via --full.
CODE_PRIORITY = [
    "CreateDensityFields.m",
    "Thompson_etal_CEE26_Figure1.m",
    "Thompson_etal_CEE26_Figure2.m",
    "Thompson_etal_CEE26_Figure3.m",
    "Thompson_etal_CEE26_Figure4.m",
    "Thompson_etal_CEE26_Figure5_6.m",
    "Thompson_etal_CEE26_SuppFigure4.m",
    "Thompson_etal_CEE26_SuppFigure6.m",
    "sshbindens_681_final_curve_090.mat",
    "possies659.mat",
    "possies681.mat",
    "w_data.mat",
    "means_SSH_1993_2017.nc",
    "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc",
    "FOCUS_data_trunc_swath_155_9.nc",
    "SWOT_Swaths_ThompsonetalCEE26.zip",
]

GLIDER_PRIORITY = [
    "SG681_SWOT_FOCUS_logfiles.tar",
    "SG659_SWOT_FOCUS_logfiles.tar",
    "SG681_SWOT_FOCUS_ncfiles.tar",
    "SG659_SWOT_FOCUS_ncfiles.tar",
]


def fetch_record(record_id: int) -> dict:
    url = f"https://zenodo.org/api/records/{record_id}"
    with urllib.request.urlopen(url, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _finalize(tmp: Path, dest: Path, expected_md5: str | None, expected_size: int | None) -> None:
    size = tmp.stat().st_size
    if expected_size is not None and size != expected_size:
        raise RuntimeError(f"size mismatch for {dest.name}: got {size}, expected {expected_size}")
    if expected_md5 is not None:
        got = md5_file(tmp)
        if got != expected_md5:
            raise RuntimeError(f"md5 mismatch for {dest.name}: got {got}, expected {expected_md5}")
    # Windows: replace can fail while AV holds the file; copy then remove.
    if dest.exists():
        dest.unlink()
    try:
        tmp.replace(dest)
    except OSError:
        shutil.copyfile(tmp, dest)
        tmp.unlink(missing_ok=True)


def _unlink_with_retry(path: Path, attempts: int = 10) -> None:
    for i in range(1, attempts + 1):
        try:
            path.unlink(missing_ok=True)
            return
        except OSError as exc:
            if i == attempts:
                raise
            print(f"[wait] cannot delete {path.name} yet ({exc}); retry {i}/{attempts}")
            time.sleep(1)


def _reset_partial_if_oversized(tmp: Path, expected_size: int | None) -> int:
    """Drop corrupt/oversized .part files so resume cannot keep appending."""
    if not tmp.exists():
        return 0
    already = tmp.stat().st_size
    if expected_size is not None and already > expected_size:
        print(f"[reset] {tmp.name} oversized ({already} > {expected_size}); deleting")
        _unlink_with_retry(tmp)
        return 0
    return already


def _acquire_download_lock(dest: Path) -> Path:
    """Prevent multiple --full runs from writing the same .part (WinError 32)."""
    lock_path = Path(str(dest) + ".lock")
    if lock_path.exists():
        try:
            old_pid = int(lock_path.read_text(encoding="utf-8").strip())
        except ValueError:
            old_pid = -1
        still_running = False
        if old_pid > 0:
            try:
                # Windows: OpenProcess via os.kill(pid, 0) is not reliable; use tasklist.
                check = subprocess.run(
                    ["tasklist", "/FI", f"PID eq {old_pid}", "/NH"],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                still_running = str(old_pid) in check.stdout
            except OSError:
                still_running = False
        if still_running:
            raise RuntimeError(
                f"another download holds {lock_path.name} (pid {old_pid}). "
                "Stop that process before re-running --full."
            )
        _unlink_with_retry(lock_path)
    lock_path.write_text(str(os.getpid()), encoding="utf-8")
    return lock_path


def _download_with_curl(url: str, tmp: Path, expected_size: int | None) -> None:
    curl = shutil.which("curl.exe") or shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl not available")
    already = _reset_partial_if_oversized(tmp, expected_size)
    if expected_size is not None and already == expected_size:
        return
    cmd = [
        curl,
        "-L",
        "--fail",
        "--retry",
        "5",
        "--retry-delay",
        "2",
        "--retry-all-errors",
        "--connect-timeout",
        "30",
        "--max-time",
        "0",
        "-o",
        str(tmp),
        url,
    ]
    if already > 0:
        o_idx = cmd.index("-o")
        cmd[o_idx:o_idx] = ["-C", "-"]
        print(f"[curl] resume from {already} bytes")
    else:
        print("[curl] fresh download")
    completed = subprocess.run(cmd, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"curl exit {completed.returncode}")


def _download_with_urllib(url: str, tmp: Path, expected_size: int | None) -> None:
    already = _reset_partial_if_oversized(tmp, expected_size)
    if expected_size is not None and already == expected_size:
        return
    headers = {"User-Agent": "swot-qg-omega-downloader/1.1"}
    if already > 0:
        headers["Range"] = f"bytes={already}-"
        print(f"[http] resume from {already} bytes")
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        status = getattr(response, "status", None) or response.getcode()
        mode = "ab" if already > 0 and status == 206 else "wb"
        if already > 0 and status == 200:
            print("[http] server ignored Range; restarting")
            already = 0
            mode = "wb"
        with tmp.open(mode) as handle:
            last_report = already
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                already += len(chunk)
                if expected_size and already > expected_size:
                    raise RuntimeError(
                        f"download exceeded expected size ({already} > {expected_size})"
                    )
                if expected_size and already - last_report >= 16 * 1024 * 1024:
                    pct = 100.0 * already / expected_size
                    print(f"[prog] {dest_name_progress(tmp)} {already}/{expected_size} ({pct:.1f}%)")
                    last_report = already


def dest_name_progress(tmp: Path) -> str:
    name = tmp.name
    return name[:-5] if name.endswith(".part") else name


def download_file(
    url: str,
    dest: Path,
    expected_md5: str | None,
    expected_size: int | None,
    max_attempts: int = 20,
) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and expected_size is not None and dest.stat().st_size == expected_size:
        if expected_md5 is None or md5_file(dest) == expected_md5:
            print(f"[skip] {dest.name}")
            return

    lock_path = _acquire_download_lock(dest)
    tmp = Path(str(dest) + ".part")
    last_error: Exception | None = None
    try:
        for attempt in range(1, max_attempts + 1):
            have = tmp.stat().st_size if tmp.exists() else 0
            print(f"[get ] {dest.name} ({expected_size} bytes; attempt {attempt}/{max_attempts}; have {have})")
            try:
                try:
                    _download_with_curl(url, tmp, expected_size)
                except Exception as curl_exc:  # noqa: BLE001
                    print(f"[warn] curl path failed ({curl_exc}); falling back to urllib")
                    _download_with_urllib(url, tmp, expected_size)
                _finalize(tmp, dest, expected_md5, expected_size)
                print(f"[ok  ] {dest.name}")
                return
            except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
                last_error = exc
                have = _reset_partial_if_oversized(tmp, expected_size)
                print(f"[retry] {dest.name}: {exc} (keeping {have} bytes)")
                time.sleep(min(2 * attempt, 30))
        raise RuntimeError(f"failed to download {dest.name} after {max_attempts} attempts") from last_error
    finally:
        _unlink_with_retry(lock_path)


def files_by_key(record: dict) -> dict[str, dict]:
    return {item["key"]: item for item in record["files"]}


def download_named(record: dict, names: list[str], dest_dir: Path, skip_large: bool, large_min_bytes: int) -> list[dict]:
    catalog = files_by_key(record)
    results = []
    for name in names:
        meta = catalog[name]
        size = int(meta["size"])
        if skip_large and size >= large_min_bytes:
            print(f"[omit] {name} ({size} bytes; pass --full to download)")
            results.append({"key": name, "status": "omitted_large", "size": size})
            continue
        checksum = meta.get("checksum", "")
        expected_md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else None
        dest = dest_dir / name
        download_file(meta["links"]["self"], dest, expected_md5, size)
        results.append({"key": name, "status": "present", "path": str(dest), "size": size, "md5": expected_md5})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="Also download large archives (SWOT zip, example w NetCDF, glider nc tarballs).",
    )
    parser.add_argument(
        "--large-min-mb",
        type=int,
        default=50,
        help="Size threshold (MB) treated as large when --full is not set.",
    )
    args = parser.parse_args()
    skip_large = not args.full
    large_min_bytes = args.large_min_mb * 1024 * 1024

    print(f"Fetching Zenodo records {ZENODO_CODE} and {ZENODO_GLIDER} ...")
    code_record = fetch_record(ZENODO_CODE)
    glider_record = fetch_record(ZENODO_GLIDER)

    code_dir = ROOT / "vendor" / "zenodo_19981347"
    glider_dir = ROOT / "data" / "glider"
    duacs_dir = ROOT / "data" / "duacs"
    derived_dir = ROOT / "data" / "derived"
    swot_dir = ROOT / "data" / "swot"

    code_results = download_named(code_record, CODE_PRIORITY, code_dir, skip_large, large_min_bytes)

    # Mirror DUACS / mean SSH / GEM / example products into data/ for a stable layout.
    mirrors = {
        "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.25deg_P1D_1731273896002.nc": duacs_dir,
        "means_SSH_1993_2017.nc": duacs_dir,
        "sshbindens_681_final_curve_090.mat": derived_dir,
        "FOCUS_data_trunc_swath_155_9.nc": derived_dir,
        "SWOT_Swaths_ThompsonetalCEE26.zip": swot_dir,
    }
    for name, target_dir in mirrors.items():
        src = code_dir / name
        if src.exists():
            target_dir.mkdir(parents=True, exist_ok=True)
            dest = target_dir / name
            if not dest.exists():
                try:
                    dest.hardlink_to(src)
                except OSError:
                    dest.write_bytes(src.read_bytes())

    glider_results = download_named(glider_record, GLIDER_PRIORITY, glider_dir, skip_large, large_min_bytes)

    manifest = {
        "paper": {
            "doi": "10.1038/s43247-026-03673-y",
            "venue": "Commun. Earth Environ. (2026)",
            "method": "SWOT local GEM + QG omega",
        },
        "sources": {
            "zenodo_code": {
                "doi": "10.5281/zenodo.19981347",
                "record_id": ZENODO_CODE,
                "title": code_record.get("title"),
            },
            "zenodo_glider": {
                "doi": "10.5281/zenodo.17252347",
                "record_id": ZENODO_GLIDER,
                "title": glider_record.get("title"),
            },
            "swot_l3": {
                "product": "SWOT Level-3 KaRIn Low Rate SSH Expert (v1.0.2)",
                "doi": "10.24400/527896/A01-2023.018",
                "provider": "AVISO/DUACS",
                "note": "FOCUS-cropped L3 Expert swaths are bundled in the code deposit zip archive.",
            },
            "duacs_l4": {
                "product": "Altimeter satellite gridded Sea Level Anomalies (DUACS)",
                "doi": "10.48670/moi-00148",
                "note": "DUACS ADT supports GEM training/validation; an example NetCDF is in the code deposit.",
            },
            "xinvert": {
                "url": "https://github.com/miniufo/xinvert",
                "doi": "10.21105/joss.05510",
                "role": "QG omega inversion",
            },
        },
        "libraries": {
            "reference_scripts": "Figure and GEM drivers shipped with Zenodo 19981347",
            "python": [
                "xinvert",
                "xarray",
                "numpy",
                "scipy",
                "numba",
                "dask",
                "netCDF4",
                "gsw",
                "glidertools",
            ],
        },
        "files": {
            "zenodo_19981347": code_results,
            "zenodo_17252347": glider_results,
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {MANIFEST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
