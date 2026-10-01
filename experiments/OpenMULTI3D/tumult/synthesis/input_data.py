"""The published input data that does not live in git, and a command line to fetch it.

Besides the atmosphere grids (:mod:`.grid_catalog`), OpenMULTI3D's larger inputs are
folders of plain files on one public Seafile share (MPCDF Keeper): the master line lists,
the public model atoms and the metal-opacity (absmet) tables. Each :class:`Dataset` names
a share folder and where its files belong under the checkout's ``input/`` directory. The
listing, planning and resumable transfer are those of :mod:`.linelist_catalog`, which also
serve the TumultStudio Preferences pages.

Command line::

    python -m tumult.synthesis.input_data list [KIND ...]
    python -m tumult.synthesis.input_data fetch KIND [KIND ...] [--file NAME ...] [--dry-run]
    python -m tumult.synthesis.input_data fetch atoms --dest DIR

``KIND`` is ``linelists``, ``atoms`` or ``absmet``. Without ``--dest`` the files go to the
dataset's folder under ``input/``; a file already complete there is not fetched again.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

from .grid_catalog import human
from . import linelist_catalog
from .linelist_catalog import ShareFolder

CHECKOUT = Path(__file__).resolve().parents[2]
_KEEPER = dict(server=linelist_catalog.SOURCE.server, token=linelist_catalog.SOURCE.token,
               citation=linelist_catalog.SOURCE.citation)


@dataclass(frozen=True)
class Dataset:
    """A share folder and its default place under ``input/``."""
    id: str
    source: ShareFolder
    folder: str            # relative to the checkout's input/ directory

    @property
    def default_dest(self):
        return CHECKOUT / "input" / self.folder


DATASETS = {d.id: d for d in (
    Dataset("linelists", linelist_catalog.SOURCE, "linelists/master_linelists"),
    Dataset("atoms", ShareFolder(id="keeper-atoms", name="Model atoms",
                                 path="/atoms/public/", **_KEEPER), "atoms/public"),
    Dataset("absmet", ShareFolder(id="keeper-absmet", name="Metal opacity tables",
                                  path="/absmet_files/", **_KEEPER), "absmet_files"),
)}


def main(argv=None):
    ap = argparse.ArgumentParser(description="fetch the published OpenMULTI3D input data")
    sub = ap.add_subparsers(dest="cmd", required=True)
    ls = sub.add_parser("list", help="list the files of one or more datasets")
    ls.add_argument("kinds", nargs="*", choices=list(DATASETS), default=[],
                    help="datasets to list; default: all")
    f = sub.add_parser("fetch", help="download datasets")
    f.add_argument("kinds", nargs="+", choices=list(DATASETS))
    f.add_argument("--dest", type=Path, help="destination folder (one dataset only)")
    f.add_argument("--file", action="append", help="a file name; default: every file")
    f.add_argument("--dry-run", action="store_true", help="print the plan only")
    args = ap.parse_args(argv)
    kinds = args.kinds or list(DATASETS)
    if args.cmd == "fetch" and args.dest is not None and len(kinds) > 1:
        ap.error("--dest takes a single dataset")

    for kind in kinds:
        dataset = DATASETS[kind]
        listing = linelist_catalog.list_files(dataset.source)
        if args.cmd == "list":
            print(f"{dataset.source.name} ({kind}):")
            for entry in listing.files:
                print(f"  {human(entry.size):>10}  {entry.name}")
            print(f"  {human(listing.size):>10}  total, {len(listing.files)} files")
            continue
        p = linelist_catalog.plan(listing, args.dest or dataset.default_dest, args.file)
        print(p.summary())
        if args.dry_run:
            continue

        def progress(done, total, name):
            sys.stdout.write(f"\r{100 * done / total:5.1f}%  {name:<40}")
            sys.stdout.flush()
        linelist_catalog.fetch(p, progress=progress)
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
