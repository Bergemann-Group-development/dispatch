"""The published model-atmosphere grids OpenMULTI3D can fetch, and how to fetch them.

Each :class:`GridSpec` names a grid, its archive and the registered-grid record that
reads it once on disk (``format`` and ``pattern``, see :mod:`.workflows`). A download is
planned before any bytes move: :func:`list_models` reads the remote grid once,
:func:`plan` takes the chosen models and the snapshot count, checks what the destination
already holds and totals the transfer, and :func:`fetch` carries the plan out, resumably,
recording the result in the grid folder's manifest. Raising the snapshot count later only
adds snapshots, since :func:`snapshot_grid.halton_order` picks nested prefixes.

The downloaded layout is what the spec's ``pattern`` scans: a flat folder of ``.mod``
files for MARCS, one folder per model holding its snapshots for the 3D grids. Only the
files the Fortran readers open are fetched (``atm3d``/``mesh`` pairs for Stagger,
``.atmos``/``.mesh`` pairs for M3DIS); a grid's EOS and opacity-binning tables are not.
One snapshot is 33 MB for Stagger and 100 to 200 MB for M3DIS, so a whole grid at one
snapshot per model is a few GB for Stagger and about 30 GB for M3DIS; a stellar
parameter box and a snapshot count keep a download small.

Command line::

    python -m tumult.synthesis.grid_catalog list
    python -m tumult.synthesis.grid_catalog models m3dis-plato
    python -m tumult.synthesis.grid_catalog plan  m3dis-plato DEST --snapshots 2 --teff 5000:6500
    python -m tumult.synthesis.grid_catalog fetch m3dis-plato DEST --model PLATO_t57.77g44.40m0.000_v1.0
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import sys

from .snapshot_grid import DEFAULT_SNAPSHOTS, halton_order
from .workflows import DEFAULT_PATTERN, M3DIS_PATTERN, STAGGER_PATTERN, compile_pattern
from . import remote

MANIFEST = "grid_manifest.json"


@dataclass(frozen=True)
class GridSpec:
    id: str
    name: str
    format: str            # the Fortran atmos_format
    pattern: str           # the registered grid's path pattern
    kind: str              # "marcs_tar" | "zenodo_zip" | "m3dis_index"
    url: str               # the archive or index
    page: str              # where the grid is published
    citation: str
    folder: str = ""       # the grid folder's name under the atmosphere root
    prefix: str = ""       # model folders to take: the zip's root folder / the index prefix

    @property
    def snapshots(self):
        return self.kind != "marcs_tar"


CATALOG = {spec.id: spec for spec in (
    GridSpec(
        id="marcs", name="MARCS (standard composition)", format="Marcs",
        pattern=DEFAULT_PATTERN, kind="marcs_tar",
        url="https://marcs.astro.uu.se/data/marcs_st_mod.tar",
        page="https://marcs.astro.uu.se/",
        citation="Gustafsson et al. 2008, A&A 486, 951", folder="MARCS"),
    GridSpec(
        id="stagger-dwarfs", name="Stagger 2.0 dwarfs", format="Stagger2",
        pattern=STAGGER_PATTERN, kind="zenodo_zip",
        url="https://zenodo.org/records/18836730/files/dwarfs.zip?download=1",
        page="https://doi.org/10.5281/zenodo.18836730",
        citation="Rodríguez Díaz et al. 2024, A&A 693, A212; Lagae et al.",
        folder="Stagger2_dwarfs", prefix="dwarfs/"),
    GridSpec(
        id="stagger-giants", name="Stagger 2.0 giants", format="Stagger2",
        pattern=STAGGER_PATTERN, kind="zenodo_zip",
        url="https://zenodo.org/records/21809779/files/giants.zip?download=1",
        page="https://doi.org/10.5281/zenodo.21809779",
        citation="Rodríguez Díaz et al. 2024, A&A 693, A212; Lagae et al.",
        folder="Stagger2_giants", prefix="giants/"),
    GridSpec(
        id="m3dis-plato", name="M3DIS PLATO", format="M3DIS",
        pattern=M3DIS_PATTERN, kind="m3dis_index",
        url="https://gemini-web.mpia.de/atmos3D/",
        page="https://gemini-web.mpia.de/atmos3D/",
        citation="Eitner et al. 2024, A&A 688, A52",
        folder="M3DIS_PLATO", prefix="PLATO_"),
)}


# --- what a grid holds ------------------------------------------------------------------

@dataclass(frozen=True)
class RemoteFile:
    """One file to fetch: ``source`` is a URL or a :class:`remote.ZipMember`; ``path`` is
    its place under the grid folder."""
    source: object
    path: str
    size: int


@dataclass
class Model:
    """One node of a remote 3D grid and the files of each of its snapshots."""
    name: str
    teff: float
    logg: float
    feh: float
    snapshots: dict = field(default_factory=dict)   # number -> [RemoteFile]

    @property
    def available(self):
        return sorted(self.snapshots)

    def picks(self, count):
        """The ``count`` snapshots a download fetches: the Halton prefix, sorted."""
        available = self.available
        return sorted(available[i] for i in halton_order(len(available))[:count])

    def size(self, count):
        """The bytes of the ``count`` snapshots :meth:`picks` chooses."""
        return sum(f.size for snap in self.picks(count) for f in self.snapshots[snap])

    def on_disk(self, dest):
        """The snapshot numbers whose files are complete under ``dest``."""
        dest = Path(dest)
        return [snap for snap, files in sorted(self.snapshots.items())
                if all(_complete(dest / f.path, f) for f in files)]


def _parse(spec, path):
    return compile_pattern(spec.pattern).search(path)


def _dedupe(models):
    """One model per node: an updated model beside its predecessor (``t45g35m30_new``)
    wins as the lexically last folder, matching the local scan's rule."""
    by_node = {}
    for model in sorted(models, key=lambda m: m.name):
        by_node[(model.teff, model.logg, model.feh)] = model
    return sorted(by_node.values(), key=lambda m: m.name)


def _zip_models(spec, archive):
    """The models of a Stagger zip from its central directory."""
    files = defaultdict(lambda: defaultdict(list))
    params = {}
    for name, member in archive.members.items():
        if not name.startswith(spec.prefix) or name.endswith("/"):
            continue
        rel = name[len(spec.prefix):]
        folder, _, filename = rel.partition("/")
        if not filename or "/" in filename:
            continue
        kind = filename.split(".", 1)[0]
        if kind not in ("atm3d", "mesh"):
            continue
        snap = _parse(spec, rel.replace("mesh.", "atm3d.", 1)) if kind == "mesh" \
            else _parse(spec, rel)
        if snap is None:
            continue
        params.setdefault(folder, snap)
        files[folder][int(snap["snap"])].append(
            RemoteFile(member, rel, member.csize))
    models = []
    for folder, snaps in files.items():
        d = params[folder]
        complete = {n: sorted(fs, key=lambda f: f.path) for n, fs in snaps.items()
                    if len(fs) == 2}
        models.append(Model(folder, d["Teff"], d["logg"], d["FeH"], complete))
    return _dedupe(models)


def _index_models(spec, cancel=None):
    """The models of the M3DIS index: folders under the prefix whose names bind the
    stellar parameters, each listed (a few at a time) for its snapshot files. ``cancel``
    is checked before every page, so a cancelled listing stops within one request."""
    models = []
    for entry in remote.list_index(spec.url):
        if not entry.is_dir or not entry.name.startswith(spec.prefix):
            continue
        d = _parse(spec, f"{entry.name}/m3dis_0.atmos")
        if d is None:
            continue
        models.append(Model(entry.name, d["Teff"], d["logg"], d["FeH"]))
    models = _dedupe(models)

    def fill(model):
        if cancel is not None and cancel.cancelled:
            return                      # the pool drains without further requests
        folder = f"{spec.url}{model.name}/"
        listing = {e.name: e for e in remote.list_index(folder) if not e.is_dir}
        for name, entry in listing.items():
            m = re.fullmatch(r"m3dis_(\d+)\.atmos", name)
            mesh = listing.get(f"m3dis_{m.group(1)}.mesh") if m else None
            if m is None or mesh is None:
                continue
            model.snapshots[int(m.group(1))] = [
                RemoteFile(folder + name, f"{model.name}/{name}", entry.size),
                RemoteFile(folder + mesh.name, f"{model.name}/{mesh.name}", mesh.size)]
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(fill, models))
    if cancel is not None:
        cancel.check()
    return models


# --- planning and fetching --------------------------------------------------------------

@dataclass
class Listing:
    """A remote grid as listed: its models, and the archive their files come from."""
    spec: GridSpec
    models: list                          # Model, sorted by name; empty for a 1D grid
    archive: object = None                # the RemoteZip the members come from
    size: int = 0                         # a 1D grid's archive size

    def find(self, name):
        return next((m for m in self.models if m.name == name), None)


def list_models(grid_id, cancel=None):
    """Read the remote grid ``grid_id`` once: its models with their available snapshots
    (three range requests for a Stagger zip, one page per model folder for M3DIS).
    ``cancel`` is a :class:`~.cancel.CancelToken` checked between requests."""
    spec = CATALOG[grid_id]
    if cancel is not None:
        cancel.check()
    if not spec.snapshots:
        return Listing(spec, [], size=remote.content_length(spec.url))
    if spec.kind == "zenodo_zip":
        archive = remote.RemoteZip(spec.url)
        if cancel is not None:
            cancel.check()
        return Listing(spec, _zip_models(spec, archive), archive)
    return Listing(spec, _index_models(spec, cancel))


@dataclass
class Plan:
    listing: Listing
    dest: Path
    snapshots: int | None                 # per model; None for a 1D grid
    models: list                          # the chosen models
    chosen: dict                          # model name -> [snapshot numbers]
    files: list                           # RemoteFile, in transfer order
    present: list                         # RemoteFile already complete on disk

    @property
    def spec(self):
        return self.listing.spec

    @property
    def total_bytes(self):
        return sum(f.size for f in self.files)

    @property
    def present_bytes(self):
        return sum(f.size for f in self.present)

    @property
    def pending_bytes(self):
        return self.total_bytes - self.present_bytes

    @property
    def pending(self):
        done = {f.path for f in self.present}
        return [f for f in self.files if f.path not in done]

    def summary(self):
        lines = [f"{self.spec.name}: {self.spec.url}", f"destination: {self.dest}"]
        if self.spec.snapshots:
            n_snaps = sum(len(v) for v in self.chosen.values())
            lines.append(f"{len(self.models)} models, {n_snaps} snapshots "
                         f"({self.snapshots} per model)")
        lines.append(f"{human(self.total_bytes)} in {len(self.files)} files, "
                     f"{human(self.present_bytes)} already present")
        return "\n".join(lines)


def human(n):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if n < 1024 or unit == "TiB":
            return f"{n:.1f} {unit}"
        n /= 1024


def _size_on_disk(f):
    return f.source.usize if isinstance(f.source, remote.ZipMember) else f.size


def _complete(path, f):
    """Whether ``path`` holds ``f``: M3DIS index sizes are rounded, so a file of about the
    right size counts; the transfer itself compares against the exact length."""
    return path.is_file() and abs(path.stat().st_size - _size_on_disk(f)) <= 0.01 * f.size


def plan(listing, dest, snapshots=None, models=None):
    """Plan fetching ``models`` (names; all of the listing by default) into ``dest``
    with ``snapshots`` per model, the count of :func:`snapshot_grid.halton_order` picks."""
    spec = listing.spec
    dest = Path(dest).expanduser()
    if not spec.snapshots:
        files = [RemoteFile(spec.url, Path(spec.url).name, listing.size)]
        done = read_manifest(dest).get("complete", False)
        return Plan(listing, dest, None, [], {}, files, files if done else [])
    count = int(snapshots if snapshots is not None else DEFAULT_SNAPSHOTS)
    if count < 1:
        raise ValueError("snapshots per model must be at least 1")
    wanted = listing.models if models is None else [
        m for m in listing.models if m.name in set(models)]
    chosen, files = {}, []
    for model in wanted:
        picks = model.picks(count)
        chosen[model.name] = picks
        for snap in picks:
            files.extend(model.snapshots[snap])
    present = [f for f in files if _complete(dest / f.path, f)]
    return Plan(listing, dest, count, wanted, chosen, files, present)


def fetch(plan, progress=None, cancel=None):
    """Carry out ``plan``: fetch every pending file, then write the manifest and return
    the registered-grid record. ``progress(done_bytes, total_bytes, label)`` is called
    as bytes arrive; ``cancel`` is a :class:`~.cancel.CancelToken`."""
    dest = plan.dest
    dest.mkdir(parents=True, exist_ok=True)
    if plan.spec.kind == "marcs_tar":
        _fetch_marcs(plan, progress, cancel)
    else:
        _fetch_snapshots(plan, progress, cancel)
    write_manifest(plan)
    return grid_record(plan.spec, dest, plan.snapshots)


def _fetch_snapshots(plan, progress, cancel):
    total = plan.total_bytes
    done = plan.present_bytes
    for f in plan.pending:
        if cancel is not None:
            cancel.check()
        label = f.path

        def report(got, _size, base=done, label=label):
            if progress is not None:
                progress(base + got, total, label)
        target = plan.dest / f.path
        if isinstance(f.source, remote.ZipMember):
            plan.listing.archive.extract(f.source, target, progress=report, cancel=cancel)
        else:
            remote.fetch_file(f.source, target, progress=report, cancel=cancel)
        done += f.size


def _fetch_marcs(plan, progress, cancel):
    """Stream the MARCS tar once, gunzipping every model into the grid folder."""
    (only,) = plan.files
    if plan.present:
        return

    def handle(member, stream):
        name = Path(member.name).name
        if not name.endswith(".mod.gz"):
            return
        target = plan.dest / name[:-3]
        if target.exists() and target.stat().st_size > 0:
            return
        with gzip.open(stream) as gz, target.open("wb") as out:
            out.write(gz.read())

    def report(got, size):
        if progress is not None:
            progress(got, size, only.path)
    remote.stream_tar(plan.spec.url, handle, progress=report, cancel=cancel)


# --- the grid folder's record -----------------------------------------------------------

def manifest_path(dest):
    return Path(dest).expanduser() / MANIFEST


def read_manifest(dest):
    path = manifest_path(dest)
    if not path.is_file():
        return {}
    return json.loads(path.read_text())


def write_manifest(plan):
    """Record the source and, per model, the snapshots fetched so far (a later fetch
    with more snapshots extends the lists)."""
    record = read_manifest(plan.dest)
    fetched = record.get("fetched", {})
    for name, picks in plan.chosen.items():
        fetched[name] = sorted(set(fetched.get(name, [])) | set(picks))
    record.update({
        "catalog_id": plan.spec.id, "name": plan.spec.name, "source": plan.spec.url,
        "format": plan.spec.format, "pattern": plan.spec.pattern,
        "updated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    })
    if plan.spec.snapshots:
        record["snapshots"] = max(plan.snapshots, int(record.get("snapshots") or 0))
        record["fetched"] = fetched
    else:
        record["complete"] = True
    manifest_path(plan.dest).write_text(json.dumps(record, indent=1, sort_keys=True))
    return record


def grid_record(spec, folder, snapshots=None):
    """The registered-grid record (see ``TumultStudio2.gridmanager``) reading ``folder``:
    the catalogue id, the spec's format and pattern, and the snapshot policy matching
    what was fetched."""
    record = {"id": spec.id, "name": spec.name, "folder": str(Path(folder).expanduser()),
              "pattern": spec.pattern, "format": spec.format}
    if spec.snapshots:
        record["snapshots"] = str(snapshots or DEFAULT_SNAPSHOTS)
    return record


def grid_folder(spec, root):
    """Where ``spec`` lives under the storage ``root``."""
    return Path(root).expanduser() / spec.folder


def _models(spec, folder):
    """The names of the models held in ``folder``, unsorted, as a lazy scan: the ``.mod``
    files of a 1D grid, the non-empty snapshot folders of a 3D one. ``os.scandir`` reads
    the directory incrementally, so a consumer that stops early pays for one read."""
    with os.scandir(folder) as entries:
        for entry in entries:
            if not spec.snapshots:
                if entry.is_file() and entry.name.endswith(".mod"):
                    yield entry.name
            elif entry.is_dir() and any(f.is_file() for f in Path(entry.path).iterdir()):
                yield entry.name


def models_on_disk(spec, root):
    """The model folders (1D: the model files) held under the storage ``root``, by name.
    A catalogue grid counts as available once this is non-empty, whether the files came
    from :func:`fetch` or elsewhere (the bundled MARCS copy has no manifest)."""
    folder = grid_folder(spec, root)
    return sorted(_models(spec, folder)) if folder.is_dir() else []


def has_models_on_disk(spec, root):
    """Whether any model of ``spec`` is under the storage ``root``. Stops at the first
    one: the registry reads this on every access, and a 1D grid holds tens of thousands
    of files."""
    folder = grid_folder(spec, root)
    return folder.is_dir() and next(_models(spec, folder), None) is not None


def delete_model(spec, root, name):
    """Remove one fetched model of ``spec`` (its folder) and drop it from the manifest."""
    folder = grid_folder(spec, root)
    target = folder / name
    if target.is_dir() and target.resolve().parent == folder.resolve():
        shutil.rmtree(target)
    record = read_manifest(folder)
    if name in record.get("fetched", {}):
        del record["fetched"][name]
        manifest_path(folder).write_text(json.dumps(record, indent=1, sort_keys=True))


def delete_grid(spec, root):
    """Remove everything fetched for ``spec`` under the storage ``root``."""
    folder = grid_folder(spec, root)
    if folder.is_dir():
        shutil.rmtree(folder)


# --- command line -----------------------------------------------------------------------

def _range(text):
    lo, _, hi = text.partition(":")
    return (float(lo) if lo else None, float(hi) if hi else None)


def _in_box(model, box):
    """Whether ``model`` lies in ``box``: ``{"teff": (lo, hi), ...}`` with either bound
    ``None`` for open; a missing key is unconstrained."""
    for key in ("teff", "logg", "feh"):
        lo, hi = box.get(key) or (None, None)
        value = getattr(model, key)
        if (lo is not None and value < lo) or (hi is not None and value > hi):
            return False
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fetch a published model-atmosphere grid.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="the grids in the catalogue")
    sub.add_parser("models", help="the models of a grid").add_argument(
        "grid", choices=sorted(CATALOG))
    for name in ("plan", "fetch"):
        p = sub.add_parser(name)
        p.add_argument("grid", choices=sorted(CATALOG))
        p.add_argument("dest", type=Path, help="the grid folder to fill")
        p.add_argument("--snapshots", type=int, default=None,
                       help=f"snapshots per model (default {DEFAULT_SNAPSHOTS})")
        p.add_argument("--model", action="append", default=[], metavar="NAME",
                       help="a model to fetch (repeatable; default: every model)")
        for key in ("teff", "logg", "feh"):
            p.add_argument(f"--{key}", type=_range, metavar="LO:HI",
                           help=f"keep models with {key} in LO:HI (either side open)")
    args = ap.parse_args(argv)
    if args.cmd == "list":
        for spec in CATALOG.values():
            print(f"{spec.id:16s} {spec.name:32s} {spec.format:9s} {spec.page}")
        return 0
    listing = list_models(args.grid)
    if args.cmd == "models":
        for model in listing.models:
            print(f"{model.name:32s} {model.teff:6.0f} {model.logg:5.2f} {model.feh:+5.2f} "
                  f"snapshots {model.available}")
        return 0
    box = {k: getattr(args, k) for k in ("teff", "logg", "feh") if getattr(args, k)}
    names = [m.name for m in listing.models
             if (not args.model or m.name in args.model) and _in_box(m, box)]
    p = plan(listing, args.dest, args.snapshots, names)
    print(p.summary())
    if args.cmd == "plan":
        for model in p.models:
            print(f"  {model.name:32s} snapshots {p.chosen[model.name]} "
                  f"of {model.available}")
        return 0
    last = [None]

    def progress(done, total, label):
        pct = 100 * done / total if total else 0
        line = f"\r{pct:5.1f}%  {human(done)}/{human(total)}  {label}"
        if line != last[0]:
            sys.stdout.write(line.ljust(100))
            sys.stdout.flush()
            last[0] = line
    record = fetch(p, progress=progress)
    print(f"\nregistered-grid record: {json.dumps(record)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
