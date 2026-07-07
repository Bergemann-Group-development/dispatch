#!/usr/bin/env python3
"""Download the OpenMULTI3D input data that is too large to live in git.

The data is hosted on MPCDF Keeper (a Seafile instance) as a public shared
folder whose layout mirrors this repository's ``input/`` directory:

    https://keeper.mpdl.mpg.de/d/3191475aebe340638422/

This script walks that share, downloads every file into ``input/`` (creating
sub-directories as needed so the contents are *merged* with the files already
tracked in git), and then recursively unpacks any ``*.zip`` it finds -- including
zips nested inside other zips.

It is safe to re-run: a file is skipped when it already exists locally with the
expected size, and a zip is skipped once it has been extracted. Only stdlib is
used, so no ``pip install`` is required.

Usage
-----
    python fetch_input_data.py              # download + unpack into ./input
    python fetch_input_data.py --dry-run    # list what would be downloaded
    python fetch_input_data.py --clean-zips # delete each zip after unpacking
    python fetch_input_data.py --dest DIR   # merge into DIR instead of ./input
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

SHARE_TOKEN = "3191475aebe340638422"
BASE_URL = "https://keeper.mpdl.mpg.de"
DIRENTS_API = f"{BASE_URL}/api/v2.1/share-links/{SHARE_TOKEN}/dirents/"
DOWNLOAD_URL = f"{BASE_URL}/d/{SHARE_TOKEN}/files/"

CHUNK = 1 << 20  # 1 MiB
RETRIES = 3


def human(nbytes: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if nbytes < 1024 or unit == "TiB":
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TiB"


def _get(url: str, params: dict) -> bytes:
    query = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    req = urllib.request.Request(f"{url}?{query}", headers={"Accept": "application/json"})
    last_err: Exception | None = None
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError) as exc:  # pragma: no cover
            last_err = exc
            if attempt < RETRIES:
                time.sleep(2 * attempt)
    raise RuntimeError(f"request failed after {RETRIES} attempts: {url} {params}: {last_err}")


def list_remote(path: str = "/") -> list[dict]:
    """Recursively list the share. Returns dicts with 'path', 'size', 'is_dir'."""
    out: list[dict] = []
    raw = _get(DIRENTS_API, {"path": path})
    for e in json.loads(raw)["dirent_list"]:
        if e["is_dir"]:
            out.extend(list_remote(e["folder_path"]))
        else:
            out.append({"path": e["file_path"], "size": int(e["size"])})
    return out


def download(remote_path: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    url = f"{DOWNLOAD_URL}?" + urllib.parse.urlencode(
        {"p": remote_path, "dl": "1"}, quote_via=urllib.parse.quote
    )
    last_err: Exception | None = None
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp, part.open("wb") as fh:
                total = int(resp.headers.get("Content-Length", 0))
                got = 0
                t0 = time.time()
                while True:
                    buf = resp.read(CHUNK)
                    if not buf:
                        break
                    fh.write(buf)
                    got += len(buf)
                    if total:
                        pct = 100 * got / total
                        sys.stdout.write(f"\r    {pct:5.1f}%  {human(got)}/{human(total)}")
                        sys.stdout.flush()
            dt = time.time() - t0
            sys.stdout.write(f"\r    done   {human(got)} in {dt:.0f}s\n")
            part.replace(dest)
            return
        except (urllib.error.URLError, TimeoutError) as exc:  # pragma: no cover
            last_err = exc
            part.unlink(missing_ok=True)
            if attempt < RETRIES:
                print(f"    retry {attempt} after error: {exc}")
                time.sleep(2 * attempt)
    raise RuntimeError(f"download failed: {remote_path}: {last_err}")


def unpack_all(root: Path, clean: bool) -> None:
    """Recursively extract every *.zip under `root` (handles zip-in-zip)."""
    while True:
        pending = [
            z for z in root.rglob("*.zip")
            if not z.with_name(z.name + ".extracted").exists()
        ]
        if not pending:
            break
        for z in pending:
            print(f"  unzip {z.relative_to(root.parent)}")
            try:
                with zipfile.ZipFile(z) as zf:
                    zf.extractall(z.parent)
            except zipfile.BadZipFile as exc:
                print(f"    WARNING: not a valid zip, skipping ({exc})")
                # mark so we don't loop forever on a corrupt file
            z.with_name(z.name + ".extracted").touch()
            if clean:
                z.unlink()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", type=Path, default=Path(__file__).resolve().parent / "input",
                    help="directory to merge the data into (default: ./input)")
    ap.add_argument("--dry-run", action="store_true", help="only list what would be downloaded")
    ap.add_argument("--clean-zips", action="store_true", help="delete each zip after unpacking")
    args = ap.parse_args()

    dest: Path = args.dest
    print(f"Listing remote share {BASE_URL}/d/{SHARE_TOKEN}/ ...")
    files = list_remote()
    total = sum(f["size"] for f in files)
    print(f"Found {len(files)} files, {human(total)} total.\n")

    if args.dry_run:
        for f in files:
            print(f"  {human(f['size']):>10}  input{f['path']}")
        return 0

    dest.mkdir(parents=True, exist_ok=True)
    downloaded = skipped = 0
    for i, f in enumerate(files, 1):
        rel = f["path"].lstrip("/")
        local = dest / rel
        marker = local.with_name(local.name + ".extracted")  # zip already unpacked+removed
        if (local.exists() and local.stat().st_size == f["size"]) or marker.exists():
            skipped += 1
            continue
        print(f"[{i}/{len(files)}] {rel}  ({human(f['size'])})")
        download(f["path"], local)
        downloaded += 1

    print(f"\nDownloaded {downloaded} file(s), skipped {skipped} already present.")
    print("Unpacking archives ...")
    unpack_all(dest, clean=args.clean_zips)
    print("Done. input/ is up to date.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
