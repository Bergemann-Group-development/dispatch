"""Everything OpenMULTI3D can download, through one interface.

Two kinds of source share it. The model-atmosphere *grids* of
:mod:`.synthesis.grid_catalog` (``marcs``, ``stagger-dwarfs``, ``stagger-giants``,
``m3dis-plato``) select their 3D models by stellar parameters; the MARCS grid is one
archive and comes whole. The *datasets* of :mod:`.synthesis.input_data` (``linelists``,
``ges``, ``atoms``, ``absmet``, ``linemasks``, ``windows``, ``kpno``, ``iag``, ``gbs``)
are folders of files on the MPCDF Keeper share, selected by file name. Each goes by
default to its folder under the checkout's ``input/`` directory, where the code and the
documentation's examples read it, and a file already complete there is not fetched
again. A source's remote listing is read once per session.

Python::

    import tumult as tm
    tm.downloads.table()                         # the catalogue
    tm.downloads.table("stagger-dwarfs")         # its models, with what is on disk
    snapshot, = tm.downloads.fetch("stagger-dwarfs", teff=5777, logg=4.4, feh=0.0)
    atom, = tm.downloads.fetch("atoms", "atom.mg86g", verbose=True)

Command line::

    python -m tumult.downloads list [ID]
    python -m tumult.downloads fetch ID [ID ...] [--file NAME ...] [--dest DIR]
    python -m tumult.downloads fetch stagger-dwarfs --teff 5777 --logg 4.4 --feh 0 --snapshots 2
"""
from __future__ import annotations

import argparse
from functools import cache
from pathlib import Path

import pandas as pd

from .synthesis import grid_catalog, input_data, linelist_catalog

GRIDS = grid_catalog.CATALOG
DATASETS = input_data.DATASETS


def _source(source_id):
    if source_id in GRIDS:
        return GRIDS[source_id]
    if source_id in DATASETS:
        return DATASETS[source_id]
    raise LookupError(f"unknown download {source_id!r}; known: {', '.join([*GRIDS, *DATASETS])}")


@cache
def _listing(source_id):
    """The remote listing of ``source_id``, read once per session."""
    if source_id in GRIDS:
        return grid_catalog.list_models(source_id)
    return linelist_catalog.list_files(DATASETS[source_id].source)


def default_dest(source_id):
    """Where ``source_id`` goes without a ``dest``: its folder under the checkout's
    ``input/`` directory."""
    if source_id in GRIDS:
        return grid_catalog.grid_folder(GRIDS[source_id], input_data.CHECKOUT / "input" / "atmos")
    return DATASETS[source_id].default_dest


def table(source_id=None, dest=None):
    """The catalogue as a DataFrame (``id``, ``kind``, ``name``, ``dest``), or the
    content of the grid or dataset ``source_id`` with what ``dest`` (default
    :func:`default_dest`) holds: for a 3D grid one row per model (``name``, ``teff``,
    ``logg``, ``feh``, ``snapshots``, ``size`` of one snapshot in bytes, ``on_disk``
    snapshots), for a dataset one row per file (``name``, ``size``, ``on_disk``), and
    for the MARCS grid its archive as one row."""
    if source_id is None:
        return pd.DataFrame(
            [(g.id, "grid", g.name, default_dest(g.id)) for g in GRIDS.values()]
            + [(d.id, "files", d.source.name, default_dest(d.id)) for d in DATASETS.values()],
            columns=["id", "kind", "name", "dest"])
    _source(source_id)
    listing = _listing(source_id)
    dest = Path(dest).expanduser() if dest is not None else default_dest(source_id)
    if source_id in DATASETS:
        return pd.DataFrame([(f.name, f.size, f.on_disk(dest)) for f in listing.files],
                            columns=["name", "size", "on_disk"])
    if not listing.spec.snapshots:
        done = grid_catalog.read_manifest(dest).get("complete", False)
        return pd.DataFrame([(Path(listing.spec.url).name, listing.size, done)],
                            columns=["name", "size", "on_disk"])
    models = listing.table()
    models["on_disk"] = [len(m.on_disk(dest)) for m in listing.models]
    return models


def fetch(source_id, *names, teff=None, logg=None, feh=None, snapshots=1, dest=None,
          verbose=False):
    """Fetch from the grid or dataset ``source_id`` into ``dest`` (default
    :func:`default_dest`) and return the local paths.

    - A 3D grid takes the model at exactly ``teff``, ``logg`` and ``feh``, the values its
      name carries (:meth:`~.grid_catalog.Listing.select` names the neighbours of a
      missing one), or every model without them, with ``snapshots`` snapshots each. The
      paths are the snapshots' atmosphere files, the ``atmos_file`` of a run.
    - The MARCS grid comes whole; the paths are its ``.mod`` files.
    - A dataset takes the files ``names``, or all of them; the paths are those files.

    Files already complete under ``dest`` are not fetched again. ``verbose`` says when
    everything is already on disk and otherwise shows a progress bar of the transfer."""
    _source(source_id)
    listing = _listing(source_id)
    dest = Path(dest).expanduser() if dest is not None else default_dest(source_id)
    star = {"teff": teff, "logg": logg, "feh": feh}
    given = [k for k, v in star.items() if v is not None]
    if source_id in DATASETS:
        if given:
            raise TypeError(f"{source_id} is a dataset of files: name them, not {', '.join(given)}")
        unknown = [n for n in names if listing.find(n) is None]
        if unknown:
            raise LookupError(f"{source_id} has no file {', '.join(unknown)}; its files: "
                           + ", ".join(f.name for f in listing.files))
        todo = linelist_catalog.plan(listing, dest, names or None)
        paths = [todo.dest / f.name for f in todo.files]
        _run(lambda progress: linelist_catalog.fetch(todo, progress=progress), todo, paths,
             listing.source.name, verbose)
        return paths
    if names:
        raise TypeError(f"{source_id} is a grid: select its models by teff, logg and feh")
    if not listing.spec.snapshots:
        if given:
            raise TypeError(f"{source_id} is fetched whole: drop {', '.join(given)}")
        todo = grid_catalog.plan(listing, dest)
        _run(lambda progress: grid_catalog.fetch(todo, progress=progress), todo, [todo.dest],
             listing.spec.name, verbose)
        return sorted(todo.dest.glob("*.mod"))
    if given and len(given) < 3:
        raise TypeError(f"{source_id} selects a model by all of teff, logg and feh")
    models = [listing.select(teff, logg, feh).name] if given else None
    todo = grid_catalog.plan(listing, dest, snapshots, models)
    paths = [todo.dest / f.path for f in todo.files
             if grid_catalog.compile_pattern(listing.spec.pattern).search(f.path)]
    _run(lambda progress: grid_catalog.fetch(todo, progress=progress), todo, paths,
         listing.spec.name, verbose)
    return paths


def _run(transfer, todo, paths, label, verbose):
    """Carry out a plan through ``transfer(progress)``, reporting as ``verbose`` asks."""
    if not todo.pending:
        if verbose:
            shown = paths if len(paths) <= 10 else [f"{len(paths)} files in {todo.dest}"]
            print("\n".join(f"already on disk: {p}" for p in shown))
        return
    if not verbose:
        transfer(None)
        return
    from tqdm.auto import tqdm      # deferred: only a verbose fetch draws a bar
    with tqdm(total=todo.total_bytes, initial=todo.present_bytes, unit="B", unit_scale=True,
              unit_divisor=1024, desc=label) as bar:
        transfer(lambda done, _total, _name: bar.update(done - bar.n))


def main(argv=None):
    ap = argparse.ArgumentParser(description="fetch OpenMULTI3D's published input data "
                                             "and model-atmosphere grids")
    sub = ap.add_subparsers(dest="cmd", required=True)
    ls = sub.add_parser("list", help="the catalogue, or the content of one grid or dataset")
    ls.add_argument("id", nargs="?", choices=[*GRIDS, *DATASETS])
    f = sub.add_parser("fetch", help="download grids or datasets")
    f.add_argument("ids", nargs="+", choices=[*GRIDS, *DATASETS])
    f.add_argument("--file", action="append", default=[], metavar="NAME",
                   help="a file of a dataset (repeatable; default: every file)")
    for key in grid_catalog.STAR_KEYS:
        f.add_argument(f"--{key}", type=float, help="the model's parameter (3D grids)")
    f.add_argument("--snapshots", type=int, default=1, help="snapshots per 3D model")
    f.add_argument("--dest", type=Path, help="destination folder (one id only)")
    args = ap.parse_args(argv)
    with pd.option_context("display.max_rows", None, "display.width", 200):
        if args.cmd == "list":
            print(table(args.id).to_string(index=False))
            return 0
    if args.dest is not None and len(args.ids) > 1:
        ap.error("--dest takes a single id")
    star = {k: getattr(args, k) for k in grid_catalog.STAR_KEYS}
    for source_id in args.ids:
        kwargs = dict(dest=args.dest, verbose=True)
        if source_id in GRIDS and GRIDS[source_id].snapshots:
            kwargs.update(star, snapshots=args.snapshots)
        fetch(source_id, *args.file, **kwargs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
