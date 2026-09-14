"""Where a line's *strength* comes from — the one number that decides which transitions
are worth naming and which head a line-list table.

The problem this solves
-----------------------
The obvious proxy, ``lam0 · 10^loggf · exp(-E_low·EK/T)`` at some fixed T, is not a strength
at all: it carries no element abundance, no ionization fraction, no molecular dissociation
and no partition function, so it cannot be compared between two species. Ranked by it, a
two-Angstrom window of the solar spectrum is headed by TiO and C2 bands that are not in the
Sun at all, and its ordering correlates with the true one at Spearman rho ~ 0.2 — no better
than chance.

HOP computes the right quantity per line while it builds the opacity::

    ll%strength(l) = max(ll%strength(l), voigt_peak*Rij / cont(f0,i,j,k))

— line-centre extinction over continuum extinction, maximised over the atmosphere. It is
dimensionless, it compares across species, and it is the criterion the RT code itself uses
to decide which lines are worth computing.

Computed, not stored
--------------------
That number used to be readable only from the ``linelist_out.txt`` a full run writes, which
meant keeping a file per reference star — up to a gigabyte for a cool one, and silently
wrong the moment the configured line list changed, since a dump only ever describes the list
its run used. ``hop_linelist_count`` / ``hop_get_linelist`` / ``hop_linelist_species`` in the
standalone C API (``HOP/standalone/hop_c_api.f90``) hand the same numbers back from a solved
state in memory. :func:`compute_band` is one EOS solve on a reference atmosphere plus one
opacity call over a band of spectrum, and it returns the strengths of every line in it.

The atmosphere is a *reference star* of :mod:`tumult.refstars` — a named set of stellar
parameters — resolved by :func:`atmosphere_for` to the nearest node of the registered MARCS
grid (:func:`set_marcs_grid`) and read at the eight depths :func:`read_marcs` samples. Eight
depths is all a strength needs, so an atmosphere is sixteen floats, read once per node and
kept. Without a MARCS grid on disk no strength can be computed and the proxy ranking below
stands in.

Measured against ``data/sed_3000-9000``'s own dump, over 5000-5002 A with that run's line
list: all 155 lines matched, **Spearman 0.978**, and the ratio to the run's values is a near
constant 1.21 (10-90% spread 1.21-1.33) — an offset from using a 1D MARCS structure where
the run had its own, which a ranking does not care about.

The GUI uses the full per-species ``master_linelists/`` folder by default. Molecular lists
dominate its size, so annotation views wider than :data:`MOLECULE_VIEW_LIMIT_AA` use the
folder's atomic files only; focused views retain every molecule. Bands are snapped to a grid
and cached process-wide, so panning inside one costs nothing at all.

Every band returns the same frame — ``ion``/``lam0``/``evi``/``loggf``/``strength``
ascending in lam0 — as :func:`tumult.linelist.load_master` does with its proxy column, so a
caller that has to fall back never branches on which of the two answered.

Threading
---------
HOP keeps process-wide global state and refuses a second live instance. The band cache is
module-level and the native call is serialised on
:data:`tumult.synthesis.runtime.HOP_LOCK`, shared by every caller in the process: two
sessions looking at the same star reuse one computation rather than fighting over the
library.
"""
from __future__ import annotations

import json
import logging
import os
import re
import threading
import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from . import refstars
from .linelist import MASTER_PATH, atomic_files, is_molecule, load_master
from tumult.synthesis.api import PROVENANCE_NAME
from tumult.synthesis.runtime import HOP_LOCK

DEFAULT_LOGG = 4.44          # filled in when a spectrum states a Teff but not these
DEFAULT_FEH = 0.0
BAND_AA = 100.0              # bands snap to this grid, so a pan usually stays inside one
BAND_STEPS = 16              # a view spans at most this many grid steps; the step doubles
                             # until it does, so wide views share a band instead of each
                             # asking for its own (see band_bounds)
MAX_BAND_AA = 6000.0         # widest band computed in one go (cost is linear in width) —
                             # wide enough that a whole optical spectrum is one band
GRID_AA = 0.005              # sampling: must resolve a line core, or peaks are under-read
DLAM_AA = 2.0                # wing search radius handed to HOP's line-list reader
MAX_BANDS = 8                # bands kept in the process-wide cache
NDEPTH = 8                   # depths of the reference atmospheres (see read_marcs)
PARAMS_NAME = "params.nml"
MOLECULE_VIEW_LIMIT_AA = 100.0
ANNOTATION_PATH = os.environ.get(
    "TUMULT_ANNOTATION_LINELIST",
    os.environ.get(
        "TUMULT_MASTER_LINELIST",
        str(Path(__file__).resolve().parent.parent
            / "input" / "linelists" / "master_linelists"),
    ),
)

_ATMOS_RE = re.compile(r"ATMOS_FILE='\s*([^']*?)\s*'", re.I)
_FEH_RE = re.compile(r"\bFEH=\s*([-+\d.Ee]+)")

# The bands are identical for every caller, so the cache is process-wide; the native call
# itself is serialised on the shared HOP_LOCK. A separate, short-lived cache lock lets a
# caller check readiness while another thread is inside HOP.
_bands_lock = threading.Lock()
_bands = {}      # (star, linelist, lo, hi) -> frame, most recent last

# the MARCS grid reference atmospheres are read from: (folder, pattern) once the GUI (or a
# script) has named it, else None; the nodes are discovered on first use and the eight
# (T, rho) pairs of a node are kept once read
_marcs_grid = None
_marcs_lock = threading.Lock()
_structures = {}  # model path -> (T, rho)


def _linelist_key(linelist):
    """Hashable identity for a file, folder, or explicit sequence of files."""
    if isinstance(linelist, (str, Path)):
        return str(linelist)
    return tuple(str(path) for path in linelist)


def usable_linelist(path):
    """Whether ``path`` is something HOP can be handed: a Tumult line-list file (its first
    line a quoted species label) or a folder holding at least one. HOP aborts the process
    on unreadable input, so a preference that is empty, missing or pointing elsewhere must
    be caught here rather than there.
    """
    if not path:
        return False
    path = Path(path).expanduser()
    if not path.exists():
        return False
    files = [path] if path.is_file() else [
        f for f in path.iterdir() if f.is_file() and not f.name.startswith(".")]
    for f in files:
        try:
            with f.open(errors="ignore") as stream:
                if stream.readline().lstrip().startswith("'"):
                    return True
        except OSError:
            continue
    return False


def annotation_linelist(linelist, lo, hi):
    """Line-list input for an annotation view spanning ``lo`` to ``hi`` Å.

    Narrow views use every file. Above :data:`MOLECULE_VIEW_LIMIT_AA`, a folder is reduced
    to its atomic files because the large molecular lists dominate construction time.
    Single-file inputs are left untouched: there is no reliable way to omit only their
    molecular blocks without producing another list.
    """
    path = Path(linelist) if isinstance(linelist, (str, Path)) else None
    if abs(float(hi) - float(lo)) > MOLECULE_VIEW_LIMIT_AA and path is not None \
            and path.is_dir():
        return atomic_files(path)
    return linelist


def _file_identity(path):
    """A cheap cache identity that changes when a local preparation input changes."""
    path = Path(path).expanduser().resolve()
    if path.is_dir():
        entries = []
        for item in sorted(p for p in path.rglob('*') if p.is_file()):
            stat = item.stat()
            entries.append((str(item.relative_to(path)), stat.st_size, stat.st_mtime_ns))
        return {'path': str(path), 'kind': 'directory', 'entries': entries}
    stat = path.stat()
    return {'path': str(path), 'kind': 'file', 'size': stat.st_size,
            'mtime_ns': stat.st_mtime_ns}


@dataclass(frozen=True)
class PreparedLineList:
    """Reusable HOP selection for one or more deterministic MARCS references."""
    folder: Path
    master: str
    references: tuple[str, ...]
    files: tuple[str, ...]
    windows: tuple[tuple[float, float], ...]
    settings: dict
    identities: tuple[dict, ...]
    aliases: dict

    def for_reference(self, reference):
        reference = str(Path(reference).expanduser().resolve())
        reference = self.aliases.get(reference, reference)
        return self.files[self.references.index(reference)]

    def to_dict(self):
        value = asdict(self)
        value['folder'] = str(self.folder)
        return value

    @classmethod
    def load(cls, path):
        path = Path(path)
        manifest = path if path.is_file() else path / 'preparation.json'
        value = json.loads(manifest.read_text())
        value['folder'] = Path(value['folder'])
        value['references'] = tuple(value['references'])
        value['files'] = tuple(value['files'])
        value['windows'] = tuple(tuple(w) for w in value['windows'])
        value['identities'] = tuple(value['identities'])
        value['aliases'] = dict(value.get('aliases') or {})
        return cls(**value)


def prepare_linelist(output_root, references, linelist, windows, *, R,
                     dlam=DLAM_AA, vmic=1.0, feh=0.0, abundances='',
                     relative_abundances=True, aliases=None, cwd=None, cancel_token=None):
    """Crop/thin once for all required MARCS references and return its artifact.

    The digest includes identities of the master list and every atmosphere, so replacing
    any input invalidates this cache without relying on a caller-selected name. HOP decides
    the retained transitions over the complete requested windows, which keeps contributors
    from other species (blends) in the same artifact. The references are pruned by
    concurrent ``hop_prune.x`` children (:func:`tumult.prune.prune_batch`).
    """
    master_id = _file_identity(linelist)
    reference_ids = tuple(_file_identity(path) for path in references)
    bounds = tuple((float(min(w)), float(max(w))) for w in windows)
    if not bounds:
        raise ValueError('Line-list preparation needs at least one window.')
    settings = {'R': float(R), 'dlam': float(dlam), 'vmic': float(vmic),
                'feh': float(feh), 'abundances': str(abundances),
                'relative_abundances': bool(relative_abundances)}
    aliases = {str(Path(key).expanduser().resolve()):
               str(Path(value).expanduser().resolve())
               for key, value in (aliases or {}).items()}
    identity = {'master': master_id, 'references': reference_ids, 'aliases': aliases,
                'windows': bounds, 'settings': settings}
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
    folder = Path(output_root).expanduser().resolve() / digest
    folder.mkdir(parents=True, exist_ok=True)
    from . import prune
    jobs = []
    files = []
    window_file = folder / 'windows.txt'
    for index, reference in enumerate(reference_ids):
        output = folder / f'reference-{index:03d}.txt'
        files.append(str(output))
        if output.is_file():
            continue
        if not jobs:
            window_file.write_text(_window_text(sorted(bounds)))
        jobs.append((output, prune.prune_namelist(
            output, Path(reference['path']).resolve(), [Path(master_id['path']).resolve()],
            window_file, R=R, dlam=dlam, vmic=vmic, feh=feh, abundances=abundances,
            relative_abundances=relative_abundances)))
    # one hop_prune.x child per uncached reference, scheduled together
    prune.prune_batch(jobs, cwd=cwd, cancel_token=cancel_token)
    artifact = PreparedLineList(folder, master_id['path'],
                                tuple(item['path'] for item in reference_ids), tuple(files),
                                bounds, settings, (master_id, *reference_ids), aliases)
    manifest = folder / 'preparation.json'
    tmp = manifest.with_name(f'{manifest.name}.{os.getpid()}.partial')
    tmp.write_text(json.dumps(artifact.to_dict(), indent=2, sort_keys=True, default=str) + '\n')
    tmp.replace(manifest)
    return artifact


def _sorted_frame(df):
    """The columns every source returns, ascending in lam0, with a usable index."""
    cols = [c for c in ("ion", "lam0", "evi", "loggf", "strength") if c in df.columns]
    out = df[cols].copy()
    out["ion"] = out.ion.astype(str).str.strip()
    return out.sort_values("lam0").reset_index(drop=True)


# --- the model atmosphere -----------------------------------------------------------------
def read_marcs(path):
    """``(T[K], rho[g/cm3])`` down a MARCS model — the two columns the EOS needs.

    The structure is in two tables keyed by the same depth index: T in the first, density in
    the second. A strength is a maximum over the cells of the state, so what the sampling has
    to preserve is *where that maximum lives* — and it lives anywhere: high-excitation lines
    peak in the deep hot layers, strong resonance and molecular ones at the very top. Only
    :data:`NDEPTH` layers are passed, spread over the model's whole index range, which is
    what makes this cheap; the cost of a band is roughly linear in the count, and the full
    56 is about 4x this for nothing. Measured over 5000-5100 Angstrom on a solar model
    against all 56 layers: Spearman 1.0000, every line matched, strengths within 10%.

    Spanning the range is the part that matters, not the count. Four layers over
    lgTau5 = -3..0 — an obvious-looking choice, and a plausible one for a *representative*
    depth — is wrong here for exactly this reason: it clips both ends of the maximum, so the
    median line comes out 2.7x too weak, the worst 30x (1000x on a 3500 K model), and a tenth
    of the lines fall under HOP's own cutoff and are not returned at all.
    """
    lines = Path(path).read_text(errors="ignore").split("\n")
    n = int(next(l for l in lines if "Number of depth points" in l).split()[0])
    i1 = next(i for i, l in enumerate(lines) if "lgTauR" in l and "Pe" in l)
    i2 = next(i for i, l in enumerate(lines) if "lgTauR" in l and "Density" in l)
    keep = np.unique(np.linspace(0, n - 1, min(NDEPTH, n)).round().astype(int))
    T = np.array([float(lines[i1 + 1 + k].split()[4]) for k in keep])
    rho = np.array([float(lines[i2 + 1 + k].split()[3]) for k in keep])
    return T, rho


def set_marcs_grid(folder, pattern=None):
    """Name the MARCS grid reference atmospheres are read from (``None`` clears it).

    ``pattern`` is the file-name pattern of :func:`tumult.synthesis.compile_pattern`, the
    MARCS one by default. The GUI calls this with its registered MARCS grid whenever the
    grid registry changes. Bands computed on an earlier grid stay cached: they are keyed by
    the model file they were computed on.
    """
    global _marcs_grid
    from tumult.synthesis.workflows import DEFAULT_PATTERN
    with _marcs_lock:
        _marcs_grid = (str(folder), pattern or DEFAULT_PATTERN) if folder else None


def marcs_ready():
    """Whether a MARCS grid has been named and its folder exists."""
    with _marcs_lock:
        grid = _marcs_grid
    return grid is not None and Path(grid[0]).expanduser().is_dir()


def _star_params(star):
    """``(teff, logg, feh, vmic)`` of a reference star, or ``None`` for an object that has
    only a key (which then identifies it)."""
    try:
        return tuple(float(getattr(star, k)) for k in refstars.PARAM_KEYS)
    except (AttributeError, TypeError, ValueError):
        return None


def _star_key(star):
    """What a band is cached under: the parameters of the star, since an editable entry
    keeps its key while its numbers change, and two entries with one set of numbers share
    the atmosphere."""
    params = _star_params(star)
    return params if params is not None else ("key", star.key)


def marcs_node(teff, logg, feh, vmic):
    """The registered MARCS model nearest to ``(teff, logg, feh, vmic)`` as a path.

    Nearest in the scaled distance of :func:`refstars.nearest` (1000 K, 1 dex, 1 dex and
    1 km/s alike), over every node the grid's pattern matches. Raises ``LookupError`` when
    no grid is named, its folder is missing or nothing in it matched.
    """
    from tumult.synthesis.workflows import discover_atmos_grid
    with _marcs_lock:
        grid = _marcs_grid
    if grid is None:
        raise LookupError("No MARCS grid is registered for reference atmospheres.")
    try:
        nodes = discover_atmos_grid(grid[0], grid[1])
    except (OSError, ValueError) as exc:
        raise LookupError(f"MARCS grid unusable for reference atmospheres: {exc}") from exc
    keys = list(nodes)
    pts = np.array([(k[0], k[1], k[2], k[3] if len(k) > 3 else vmic) for k in keys], float)
    target = np.array([teff, logg, feh, vmic], float)
    scale = np.array([1000.0, 1.0, 1.0, 1.0])
    d = (((pts - target) / scale) ** 2).sum(axis=1)
    path = nodes[keys[int(np.argmin(d))]]
    return Path(path[0] if isinstance(path, (list, tuple)) else path)


def atmosphere_for(star):
    """``(model path, T, rho)`` for a reference star: its nearest MARCS node, read once."""
    params = _star_params(star)
    if params is None:
        raise LookupError(f"{star!r} carries no stellar parameters.")
    path = marcs_node(*params)
    with _marcs_lock:
        hit = _structures.get(path)
    if hit is None:
        hit = read_marcs(path)
        with _marcs_lock:
            _structures[path] = hit
    return path, hit[0], hit[1]


# --- the computation ----------------------------------------------------------------------
def band_bounds(lo, hi):
    """``(lo, hi)`` of the cache band covering ``[lo, hi]``, snapped to :data:`BAND_AA`.

    Snapping is what makes panning free: neighbouring views ask for the same band and hit the
    cache. The step is :data:`BAND_AA` while the view spans at most :data:`BAND_STEPS` of
    them, and doubles until it does — so a narrow view keeps the fine grid it always had,
    while every whole-spectrum view lands on the same coarse band and computes once.

    A centred crop to :data:`MAX_BAND_AA` used to apply from 1000 Å up, on the reasoning that
    a frame holds only a few tens of names anyway. That is true, but they belong *spread
    across the view*: cropping to the middle left the annotation overlay drawing every name
    it had in a narrow clump at the centre of a full-spectrum plot, with the rest of the
    range bare. The cap is now wide enough to cover an optical spectrum end to end and only
    bites on something far wider.
    """
    lo, hi = float(min(lo, hi)), float(max(lo, hi))
    if hi - lo > MAX_BAND_AA:
        mid = 0.5 * (lo + hi)
        lo, hi = mid - MAX_BAND_AA / 2, mid + MAX_BAND_AA / 2
    step = BAND_AA
    while hi - lo > BAND_STEPS * step:
        step *= 2
    return np.floor(lo / step) * step, np.ceil(hi / step) * step


def _hop_band(T, rho, linelist, lo, hi):
    """HOP line strengths in ``[lo, hi]`` Å on one structure; the caller holds HOP_LOCK."""
    import tumult as tm
    from tumult.physics import air2vac
    lam = air2vac(np.arange(lo, hi, GRID_AA))
    hop = tm.HopEOS(feh=0.0)
    try:
        hop.load_state(T, rho=rho)
        df = hop.line_strengths(lam, linelist, dlam=DLAM_AA, vmic=1.0)
    finally:
        hop.close()
    return _sorted_frame(df[(df.lam0 >= lo) & (df.lam0 <= hi)])


def band_on_atmosphere(atmosphere, linelist, lo, hi):
    """Strengths of every line in ``[lo, hi]`` Å from HOP on the MARCS model ``atmosphere``.

    The uncached counterpart of :func:`compute_band` for a caller that names the model
    itself rather than a reference star resolved through the registered grid.
    """
    T, rho = read_marcs(atmosphere)
    with HOP_LOCK:
        return _hop_band(T, rho, linelist, *band_bounds(lo, hi))


def compute_band(star, linelist, lo, hi):
    """Strengths of every line in ``[lo, hi]`` Å, from HOP on the ``star`` atmosphere.

    ``star`` is a :class:`tumult.refstars.RefStar`; its atmosphere is the MARCS node
    :func:`atmosphere_for` resolves it to. Blocking, and seconds at worst for a wide band,
    so call it off the document thread. Serialised on
    :data:`tumult.synthesis.runtime.HOP_LOCK` because HOP is a process-wide singleton. A
    separate short-lived cache lock lets the UI check readiness while another thread
    computes a band. ``None`` is cached for a band HOP could not compute, and for a star no
    MARCS grid can resolve.
    """
    key = (_star_key(star), _linelist_key(linelist), lo, hi)
    with _bands_lock:
        if key in _bands:
            return _bands[key]
    with HOP_LOCK:
        with _bands_lock:
            if key in _bands:
                return _bands[key]
        try:
            _model, T, rho = atmosphere_for(star)
            frame = _hop_band(T, rho, linelist, lo, hi)
        except LookupError as exc:
            logging.warning("Line strengths: %s", exc)
            frame = None
        except Exception:
            logging.exception("Line strengths: HOP failed for %s over %.1f-%.1f",
                              star.key, lo, hi)
            frame = None
        with _bands_lock:
            _bands[key] = frame
            while len(_bands) > MAX_BANDS:
                _bands.pop(next(iter(_bands)))
        return frame


def band_ready(star, linelist, lo, hi):
    """Whether :func:`compute_band` would return without computing anything."""
    with _bands_lock:
        return (_star_key(star), _linelist_key(linelist), lo, hi) in _bands


def cached_band(star, linelist, lo, hi):
    """The computed band if the cache holds it, else None, never computing.

    What a repaint needs: readiness and the frame in one step, since between a
    :func:`band_ready` check and a separate read another thread can evict the band.
    """
    with _bands_lock:
        return _bands.get((_star_key(star), _linelist_key(linelist), lo, hi))


# --- cropping a line list to an atmosphere ------------------------------------------------
def crop_linelist(out_path, atmos_file, linelist, windows, *, R, dlam=DLAM_AA, vmic=1.0,
                  feh=0.0, abundances="", relative_abundances=True, base=None, cwd=None,
                  cancel_token=None):
    """Write the lines a run on ``atmos_file`` would actually use, as an input line list.

    The selection a run makes is a property of its atmosphere, its windows and the
    abundances it synthesises at -- not of anything else in the namelist. So it can be made
    once per atmosphere node and handed to every run on that node as its ``linelist_file``,
    which is what this produces. The runs then read a few hundred lines instead of the whole
    master list, and none of them depends on another's output.

    The selection is made by ``hop_prune.x`` in a child process (:mod:`tumult.prune`), from
    the same namelist groups the runs read. Two crops are applied, both HOP's own.
    ``windows`` (air Angstrom, expanded by ``dlam`` and merged) bounds which lines are read
    at all; then the opacity over the windows' constant-``R`` grid marks the lines that
    clear HOP's negligibility test on this structure, and only those are written. ``R`` and
    ``dlam`` therefore have to match what the runs will use, or the crop drops contributors
    they would have had. ``base`` carries the runs' composition, background and line-list
    settings into the selection; without it the executable's defaults apply, which are the
    run defaults.

    Crop at the *strongest* abundance a caller will synthesise (``abundances``, in the
    ``Fe=+0.40`` form, with ``relative_abundances`` saying which convention): the element's
    own lines grow with its abundance, so the top rung's selection covers every weaker rung
    for them. Lines of other species can move the other way when the override changes the
    electron density, which is what keeps the run's own selection from the pruned list in
    place as the final word.

    Only :data:`tumult.prune.SAMPLE_COUNT` depths of the model are evaluated, which is what
    makes this cheap next to the runs it feeds -- the strength is a maximum over depth, so
    what the sampling has to do is span the structure, not resolve it (see
    :func:`read_marcs`).
    """
    from . import prune

    bounds = [(float(min(w)), float(max(w))) for w in windows]
    if not bounds:
        raise ValueError("crop_linelist needs at least one window")
    bounds.sort()
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    window_file = out_path.with_suffix('.windows.txt')
    window_file.write_text(_window_text(bounds))
    nml = prune.prune_namelist(out_path, Path(atmos_file).resolve(), [Path(linelist).resolve()],
                               window_file, R=R, dlam=dlam, vmic=vmic, feh=feh,
                               abundances=abundances, relative_abundances=relative_abundances,
                               base=base)
    prune.prune_batch([(out_path, nml)], cwd=cwd, cancel_token=cancel_token)
    return out_path


def _window_text(bounds):
    """The ``window_file`` body, as ``spectrum_params`` reads it (air Angstrom)."""
    return "\n".join(["# blue_AA  red_AA", *(f"{lo:.4f}  {hi:.4f}" for lo, hi in bounds)]) + "\n"


def run_params(folder):
    """``(teff, logg, feh)`` of a run folder, for picking the reference atmosphere.

    A run the synthesis API launched states them in its own ``provenance.json``. Any other
    run folder -- one produced by hand, outside the cache -- has only its ``params.nml``,
    whose atmosphere file name carries the real values; the folder name is rounded (``t58``
    for 5750 K) and is not trusted.
    """
    try:
        record = json.loads((Path(folder) / PROVENANCE_NAME).read_text())
        star = record["star"]
        return (float(star["teff"]), float(star["logg"]), float(star["feh"]))
    except (OSError, ValueError, KeyError, TypeError):
        pass
    try:
        txt = (Path(folder) / PARAMS_NAME).read_text(errors="ignore")
    except OSError:
        return None
    m = _ATMOS_RE.search(txt)
    if not m:
        return None
    a = re.search(r"[ps](\d+)_g([+-]?[\d.]+)", Path(m.group(1)).name)
    if not a:
        return None
    feh = _FEH_RE.search(txt)
    return (float(a.group(1)), float(a.group(2)),
            float(feh.group(1)) if feh else DEFAULT_FEH)


# --- source precedence: which source answers for an object and a window -------------------
# Two sources, both handing back the same frame so callers never branch on which answered:
#
# 1. **Computed** — HOP on a reference atmosphere, with the configured line list. Which
#    atmosphere is the caller's choice: a named class, or ``auto``, which follows the
#    spectrum's own (Teff, log g, [Fe/H]) — from a synthesis's parameters or a run folder's
#    namelist — to the nearest class, and falls back to the solar one for a spectrum that
#    says nothing.
# 2. **The master line list and its proxy** — when HOP cannot be reached (no ``libhop``) or
#    the band has not been computed yet. The only source whose ordering is not physical, and
#    it is labelled ``approximate`` wherever it is shown.

class Source:
    """One source's answer: the lines, what to call it, and whether it can be trusted to rank.

    ``key`` is what the frame *is*, as opposed to ``label``, which is what a footer calls it.
    A consumer caching the frame has to key on this: the label names the star but not the
    band, so two answers a pan apart look identical by it while holding different lines.
    """

    __slots__ = ("lines", "label", "kind", "params", "key")

    def __init__(self, lines, label, kind, params=None, key=None):
        self.lines = lines
        self.label = label
        self.kind = kind       # "computed" | "proxy"
        self.params = params
        self.key = key if key is not None else (kind, label)

    @property
    def approximate(self):
        return self.kind == "proxy"

    def __repr__(self):
        n = 0 if self.lines is None else len(self.lines)
        return f"<Source {self.kind} {self.label!r} n={n}>"


class LineStrengths:
    """The strength sources in precedence order, shared by the annotations and the
    line-list panel."""

    def __init__(self, master_path=MASTER_PATH, linelist=None,
                 refstar=refstars.AUTO_KEY):
        self.master_path = str(master_path)
        # HOP's input may be the complete per-species folder while master_path remains the
        # compact fallback/search catalogue used until a computed band is ready.
        self.linelist = str(linelist or master_path)
        self.linelist_usable = usable_linelist(self.linelist)
        # the chosen reference atmosphere, or AUTO_KEY to follow each spectrum's own. Via
        # the setter, because this normally arrives from a stored preference and a key that
        # has since left the table has to land on auto, not on an unrelated star.
        self.refstar = refstars.AUTO_KEY
        self.set_refstar(refstar)
        self._proxy = None
        self._atomic_proxy = None

    # --- configuration ---------------------------------------------------------------
    def set_master_path(self, path):
        self.master_path = str(path)
        self.linelist = str(path)
        self.linelist_usable = usable_linelist(self.linelist)
        self._proxy = None
        self._atomic_proxy = None

    def set_refstar(self, key):
        """Choose the reference atmosphere, by :mod:`refstars` key or ``auto``.

        Returns what stuck, so a caller can persist exactly what it set. An unknown key
        falls back to ``auto`` rather than to some arbitrary star: a stale preference should
        follow the spectrum, not silently rank it as something it is not.
        """
        key = str(key)
        self.refstar = key if (key == refstars.AUTO_KEY or key in refstars.by_key()) \
            else refstars.AUTO_KEY
        return self.refstar

    # --- what star is this? ------------------------------------------------------------
    def params_of(self, obj):
        """``(teff, logg, feh)`` stated by ``obj`` — its synthesis parameters or its run's
        namelist — or ``None`` when it says nothing about its star."""
        meta = getattr(obj, "meta", None) or {}
        stellar = meta.get("stellar") if isinstance(meta, dict) else None
        if isinstance(stellar, dict) and stellar.get("teff"):
            try:
                return (float(stellar["teff"]),
                        float(stellar.get("logg", DEFAULT_LOGG)),
                        float(stellar.get("feh", DEFAULT_FEH)))
            except (TypeError, ValueError):
                pass
        source = getattr(obj, "source", None)
        if source and len(source) > 1 and source[0] == "run" and source[1]:
            found = run_params(source[1])
            if found:
                return found
        return None

    def star_for(self, obj):
        """The :class:`refstars.RefStar` the strengths for ``obj`` are computed on.

        A chosen class wins outright — that is what choosing one is for, and it lets a
        spectrum be ranked deliberately against something other than itself. Under ``auto``
        the spectrum's own parameters pick the nearest class, and a spectrum that states none
        gets the solar one, which is the honest default for an unlabelled observation. An
        entry that has left the catalogue since it was chosen resolves as ``auto`` does.
        ``None`` only on an empty catalogue.
        """
        if self.refstar != refstars.AUTO_KEY and self.refstar in refstars.by_key():
            return refstars.get(self.refstar)
        params = self.params_of(obj)
        return refstars.nearest(*params) if params else refstars.get(refstars.DEFAULT_KEY)

    def request_for(self, obj, lo, hi):
        """Resolved ``(star, line-list input, cache band)`` for an annotation view."""
        return (self.star_for(obj), annotation_linelist(self.linelist, lo, hi),
                band_bounds(lo, hi))

    # --- source precedence -----------------------------------------------------------
    def lines_for(self, obj, lo, hi, allow_compute=False):
        """The best :class:`Source` for ``obj`` over ``[lo, hi]`` Å.

        With ``allow_compute`` false — what a repaint on the document thread passes — a band
        that has not been computed is not waited for; the proxy answers instead and
        :meth:`warm` does the computing in the background, after which the next repaint gets
        the real thing.
        """
        if not self.linelist_usable:
            return self._proxy_source()
        star, linelist, band = self.request_for(obj, lo, hi)
        atomic_only = _linelist_key(linelist) != _linelist_key(self.linelist)
        if not self.computable:
            return self._proxy_source(atomic_only=atomic_only)
        if allow_compute:
            frame = compute_band(star, linelist, *band)
        else:
            # Never enter compute_band from a repaint, even if another thread
            # evicts a previously ready band between the readiness check and read.
            frame = cached_band(star, linelist, *band)
        if frame is not None and len(frame):
            scope = " · atomic lines" if atomic_only else ""
            return Source(frame, f"HOP · {star.label} ({star.summary}){scope}",
                          "computed", (star.teff, star.logg, star.feh),
                          key=("hop", _star_key(star), _linelist_key(linelist), band))
        return self._proxy_source(atomic_only=atomic_only)

    def warm(self, obj, lo, hi):
        """Compute what :meth:`lines_for` wants, blocking. Call off the document thread."""
        return self.lines_for(obj, lo, hi, allow_compute=True)

    def pending(self, obj, lo, hi):
        """True when :meth:`warm` would change the answer — i.e. the band is not computed."""
        if not self.linelist_usable or not self.computable:
            return False
        star, linelist, band = self.request_for(obj, lo, hi)
        return not band_ready(star, linelist, *band)

    @property
    def computable(self):
        """Whether a band can be computed at all: a MARCS grid to read an atmosphere from
        and at least one reference star to read it for."""
        return marcs_ready() and bool(refstars.stars())

    def describe(self, obj, lo, hi):
        """One short line naming the source in play, for the panel footer."""
        src = self.lines_for(obj, lo, hi)
        if src.approximate and self.linelist_usable and not self.computable:
            why = ("no reference stars" if not refstars.stars()
                   else "no MARCS grid for the reference atmosphere")
            return f"strengths: approximate ({why})"
        return ("strengths: approximate (computing…)" if src.approximate
                else f"strengths: {src.label}")

    # --- the proxy source ----------------------------------------------------------------
    def _proxy_source(self, atomic_only=False):
        """The configured master list, with whatever strength column it has.

        Normally the plain 6-column form, so this is the Boltzmann proxy ``load_master``
        falls back to — the one ordering that does not survive a change of species. A list
        that happens to carry real strengths (a ``linelist_out.txt``) keeps them, and then
        the source stops being approximate.
        """
        if self._proxy is None:
            try:
                master = load_master(self.master_path)
            except Exception:
                logging.exception("Line strengths: could not read %s", self.master_path)
                self._proxy = Source(None, "unavailable", "proxy")
            else:
                proxy = bool(master.attrs.get("strength_is_proxy", True))
                name = Path(self.master_path).name
                self._proxy = Source(_sorted_frame(master),
                                     f"{name} (approximate)" if proxy else name,
                                     "proxy" if proxy else "computed",
                                     key=("master", self.master_path))
        if not atomic_only or self._proxy.lines is None:
            return self._proxy
        if self._atomic_proxy is None:
            keep = ~self._proxy.lines.ion.map(is_molecule)
            self._atomic_proxy = Source(
                self._proxy.lines[keep].reset_index(drop=True),
                f"{self._proxy.label} · atomic lines", self._proxy.kind,
                self._proxy.params, key=(*self._proxy.key, "atomic"))
        return self._atomic_proxy
