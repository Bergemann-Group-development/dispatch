"""Public synthesis entry points (signatures unchanged from the pre-split monolith):
``build_stellar_grid``, ``compute_spectrum``, ``adaptive_fit``,
``get_NLTE_correction(s)``, and the bundled namelist ``preset`` lookup.

Orchestration only — atmosphere/vmic node selection lives in
:mod:`tumult.synthesis.atmos_select`, parallel job execution in
:mod:`tumult.synthesis.executor`, cancellation in :mod:`tumult.synthesis.cancel`.
"""
from collections import defaultdict
import copy
import hashlib
import itertools
import logging
import json
import re
import os
import numpy as np
import pandas as pd
from pathlib import Path
from tumult.auxfuncs import update_dict
from tumult.decorators import instance_cache
import tumult as tm
from tumult.flux import write_set_manifest
from tumult.launch import prepare_input

from .atmos_select import select_nodes
from .executor import _make_da_jobs, _run_batch, job_output

# Bundled namelist presets (base / nlte / departures / corr_ref), co-located with the core so
# both the library and any front-end (e.g. TumultStudio) reference one canonical location.
PRESETS_DIR = Path(__file__).resolve().parents[1] / "presets"


def preset(name):
    """Absolute path to a bundled namelist preset, e.g. ``preset("base.nml")``."""
    return PRESETS_DIR / name


def _is_nlte(settings):
    """True when ``settings`` describe an NLTE solve (an atom plus iterations)."""
    if 'atom_params' not in (settings or {}):
        return False
    maxiter = (settings.get('rt_iteration_params', {}) or {}).get(
        'maxiter', (settings.get('m3d_params', {}) or {}).get('maxiter', 0))
    return maxiter > 0


# --- shared building blocks (used by both build_stellar_grid and compute_spectra) --------
# These are faithful extractions of the inline logic they replace; the serial path stays
# bit-identical (guarded by the golden tests), and the batch path reuses them so a batched
# variant issues exactly the runs its serial compute_spectrum counterpart would.

def _run_name(teff_i, logg_i, feh_i, v):
    return f't{teff_i:04.0f}g{logg_i*10:02.0f}m{feh_i*10:+03.0f}v{round(v*10):02.0f}'


def _node_name(jkey, teff_i, logg_i, feh_i, v):
    """Run directory name of one node: the 1D name keeps its historical form (and cache
    identity); a snapshot node (``v is None``) names its snapshot instead of a vmic level
    and keeps two log g decimals, since 3D grids sit at quarter steps."""
    if v is None:
        # no dots: the launcher derives the log name and output folder by suffix
        tag = re.sub(r'[^A-Za-z0-9_+-]', '_', str(jkey[1]))
        return f't{teff_i:04.0f}g{logg_i*100:03.0f}m{feh_i*10:+03.0f}_{tag}'
    return _run_name(teff_i, logg_i, feh_i, v)


def _node_atmos(atmos_path, feh_i, v):
    """The per-node ``atmos_params`` override. A snapshot's velocity field replaces
    microturbulence, so no ``vmic`` is written for it."""
    override = {'atmos_file': atmos_path, 'feh': feh_i}
    if v is not None:
        override['vmic'] = v
    return override


def _dep_jobs(base, nodes, departures_top):
    """Departure-coefficient jobs (one per atmosphere node) for the NLTE two-step path: no
    line list, no output spectrum (``rt_output_params angle_quad='none'``), saved via the
    'save' format."""
    jobs = []
    for jkey, teff_i, logg_i, feh_i, v, _w, atmos_path, absmet_file in nodes:
        dep_nml = copy.deepcopy(base)
        dep_nml.pop('linelist_params', None)
        dep_nml.pop('spectrum_params', None)
        dep_nml = update_dict(dep_nml, {
            'atmos_params': _node_atmos(atmos_path, feh_i, v),
            'opac_table_params': {'absmet_file': absmet_file},
            'rt_iteration_params': {'maxiter': 99},
            'rt_output_params': {'maxiter': 3, 'angle_quad': 'none'},
        })
        jobs.append((jkey, Path(departures_top) / _node_name(jkey, teff_i, logg_i, feh_i, v),
                     dep_nml))
    return jobs


def _spec_jobs(base, nodes, top, nlte, dep_dirs, node_linelists=None):
    """Spectrum jobs (one per node): LTE reads the marcs atmosphere directly; NLTE reads the
    saved atmosphere + precomputed departures ('save' format, maxiter=0). Returns
    ``(jobs, job_keys, job_weights, job_params)``.

    ``node_linelists`` (optional) maps a node key to the line-list files that node's runs read
    instead of the shared one -- the per-node pruned list from :func:`_node_linelists`. It
    replaces ``linelist_params.linelist_file``, so whatever produced the mapping owns the cache
    identity of these runs (see :func:`_prune_digest`)."""
    jobs, job_keys, job_weights, job_params = [], [], [], []
    for jkey, teff_i, logg_i, feh_i, v, weight, atmos_path, absmet_file in nodes:
        if nlte:
            atmos_override = dict(_node_atmos(dep_dirs[jkey], feh_i, v),
                                  atmos_format='save', read_departures=True)
            extra = {'rt_iteration_params': {'maxiter': 0}}
        else:
            atmos_override = _node_atmos(atmos_path, feh_i, v)
            extra = {}
        if node_linelists:
            extra = dict(extra,
                         linelist_params={'linelist_file': list(node_linelists[jkey])})
        nml_i = update_dict(copy.deepcopy(base), {
            'atmos_params': atmos_override,
            'opac_table_params': {'absmet_file': absmet_file},
            **extra,
        })
        jobs.append((jkey, Path(top) / _node_name(jkey, teff_i, logg_i, feh_i, v), nml_i))
        job_keys.append(jkey)
        job_weights.append(weight)
        job_params.append((teff_i, logg_i, feh_i, 0.0 if v is None else v))
    return jobs, job_keys, job_weights, job_params


def _prune_digest(base, windows, prune):
    """Cache identity of a per-node pruned line list.

    Everything the pruning depends on except the atmosphere node, which the file name carries:
    the master list, the windows, the skirt, the sampling the selection was made at and the
    abundances it was made for. This goes into the *runs'* digest as well as into the name of
    the pruned lists' directory, because two runs on the same node with the same namelist but
    different pruned lists are different runs -- without it the second would be served the
    first's cache.
    """
    files = (base.get('linelist_params') or {}).get('linelist_file') or []
    R = (base.get('spectrum_params') or {}).get('R_sampling')
    return hashlib.sha1(json.dumps(
        [[str(f) for f in files], [[float(min(w)), float(max(w))] for w in windows], prune,
         None if R is None else float(R)],     # one identity whether R came as int or float
        sort_keys=True, default=str).encode()).hexdigest()[:8]


# composition_patch_mod's default line-list read margin (&linelist_params dlam)
DEFAULT_DLAM_AA = 100.0


def _resolved_masters(files, launch_kwargs):
    """The configured line-list files as absolute paths, in the order HOP loads them, so
    a transition listed twice resolves the same way it would for a run."""
    cwd = Path((launch_kwargs or {}).get('cwd') or '.')
    return [(cwd / f).resolve() for f in files]


def _node_linelists(top, base, nodes, windows, prune, launch_kwargs, preparation=None,
                    cancel_token=None):
    """Prune the line list once per atmosphere node, and map each node key to its own file.

    A run's line selection depends on its atmosphere, its windows and its abundances, so a
    ladder of abundance rungs on one node re-reads and re-selects the whole master list at
    every rung. Pruning once per node, at the composition the ladder is centred on, makes
    the selection a property of the node instead, and leaves the runs independent of each
    other: no run waits for another's ``linelist_out.txt``. A rung some way off that
    composition sees the selection shift slightly with it, which is the approximation one
    selection per node buys.

    The selection is made by one ``hop_prune.x`` child per uncached node
    (:mod:`tumult.prune`), fed the runs' own composition, background and line-list settings
    from ``base`` and scheduled together on the workflow's CPU budget; ``cancel_token``
    terminates live children and none of them publishes a partial file. The children read
    the master lists themselves: HOP bisects each block for the windows it was given
    (``HOP/linelist_reader_mod.f90``), so finding the records costs the selection rather
    than the list.

    ``prune`` is the caller's descriptor: ``skirt`` in Angstrom (how far outside a window a
    line may still contribute; 0 or absent disables pruning altogether) and the
    ``abundances`` / ``relative_abundances`` the selection is made at, which is the
    composition itself -- a ladder's rungs share their node's selection rather than each
    getting its own. Hydrogen lines, whose Stark wings reach far beyond the skirt, bypass
    the prune and are read as the master lists give them. Each node then maps to the list
    of its files.

    Without pruning there is nothing to select per node: every run reads the shared list
    and keeps the records its own windows need, and this returns ``None``.

    A prune costs a fraction of a second against the several seconds of the runs it feeds,
    and is cached on disk under the digest, so a node recurring in a later iteration or
    session pays nothing.
    """
    if preparation is not None:
        mapping = {}
        for jkey, _teff, _logg, _feh, _v, _weight, atmos_path, _absmet in nodes:
            try:
                mapping[jkey] = [preparation.for_reference(atmos_path)]
            except ValueError as exc:
                raise ValueError(f'Prepared line list has no entry for {atmos_path}') from exc
        return mapping
    linelist = base.get('linelist_params') or {}
    files = [f for f in linelist.get('linelist_file') or [] if f]
    if not files:
        return None
    if not prune or not float(prune.get('skirt') or 0.0) > 0.0:
        return None                   # nothing to select per node; every run shares the list
    R = float((base.get('spectrum_params') or {}).get('R_sampling') or 0.0)
    if R <= 0:
        return None
    from .. import prune as hop_prune
    from ..linestrength import _window_text

    skirt = float(prune['skirt'])
    padded = sorted((float(min(w)) - skirt, float(max(w)) + skirt) for w in windows)
    cwd = Path((launch_kwargs or {}).get('cwd') or '.')
    out_dir = (cwd / top / '_linelists' / _prune_digest(base, windows, prune)).resolve()
    # HOP loads the files in this order for the prune as a run would for itself, so the
    # written list resolves a transition listed twice the same way; absolute paths
    # replace cwd
    masters = _resolved_masters(files, launch_kwargs)
    out_dir.mkdir(parents=True, exist_ok=True)
    window_file = out_dir / 'windows.txt'
    relative = bool(prune.get('relative_abundances', True))
    abundances = str(prune.get('abundances') or '')
    # hydrogen is read unpruned beside the selection, so the selection leaves it out
    hydrogen = hop_prune.hydrogen_lists(masters, out_dir)
    jobs, mapping = [], {}
    for jkey, teff_i, logg_i, feh_i, v, _w, atmos_path, _absmet in nodes:
        if v is None:
            # hop_prune.x samples 1D structures; a snapshot node needs the request path's
            # MARCS preparation references (see workflows.resolve_request) rather than a guess.
            raise ValueError('Line-list pruning on a 3D snapshot grid needs a MARCS '
                             'preparation grid; pass a prepared line list or disable pruning.')
        atmos = (atmos_path if Path(atmos_path).is_absolute() else cwd / atmos_path).resolve()
        path = out_dir / f'{_run_name(teff_i, logg_i, feh_i, v)}.txt'
        mapping[jkey] = [*hydrogen, str(path)]
        if path.is_file():
            continue
        if not jobs:
            window_file.write_text(_window_text(padded))
        jobs.append((path, hop_prune.prune_namelist(
            path, atmos, masters, window_file, R=R, dlam=skirt, vmic=v, feh=feh_i,
            abundances=abundances, relative_abundances=relative, base=base,
            exclude_elements='H')))
    hop_prune.prune_batch(jobs, cwd=cwd, cancel_token=cancel_token)
    return mapping


def _is_snapshot_key(key):
    """Snapshot-grid nodes are keyed ``((teff, logg, feh), snapshot_tag)``; 1D nodes carry
    the vmic level in that slot (see :func:`select_nodes`)."""
    return isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], str)


class _SnapshotMeanView:
    """One snapshot-grid node as :class:`SpectrumGrid` expects it: the equal-weight mean
    over the node's snapshot runs, so the grid sees one spectrum per stellar point."""
    def __init__(self, runs):
        self.runs = list(runs)
        self.has_spectrum = all(getattr(run, 'has_spectrum', False) for run in self.runs)
        self.mode = getattr(self.runs[0], 'mode', 'LTE')
        self.folder = self.runs[0].folder

    def crop(self, *args, **kwargs):
        return tm.Spectrum.combine([run.crop(*args, **kwargs) for run in self.runs])


def _collapse_snapshots(runs, job_keys, job_weights, job_params):
    """Merge the per-snapshot nodes of a snapshot grid onto their stellar point. Snapshot
    node keys carry the snapshot tag as their second element; their weights (the spatial
    weight shared equally in :func:`select_nodes`) add back up to the node's."""
    spectra, weights, params = [], [], []
    grouped = {}    # spatial key -> index into the lists
    for key, weight, point in zip(job_keys, job_weights, job_params):
        if not _is_snapshot_key(key):
            spectra.append(runs[key])
            weights.append(weight)
            params.append(point)
            continue
        index = grouped.get(key[0])
        if index is None:
            grouped[key[0]] = len(spectra)
            spectra.append(_SnapshotMeanView([runs[key]]))
            weights.append(weight)
            params.append(point)
        else:
            spectra[index].runs.append(runs[key])
            weights[index] += weight
    return spectra, weights, params


def _assemble_grid(runs, star, vmic, job_keys, job_weights, job_params, vmic_points):
    """Load the finished node runs into a tm.grid at ``star`` (separable vmic axis when
    ``vmic_points > 2``, else linear bracket with stored interpolation weights)."""
    teff, logg, feh = star[:3]
    # snapshot nodes (keyed by snapshot tag) have no vmic lattice to collapse
    sep_vmic = vmic_points > 2 and not any(_is_snapshot_key(k) for k in job_keys)
    spectra, weights, params = _collapse_snapshots(runs, job_keys, job_weights, job_params)
    param_dict = {
        'teff': [p[0] for p in params],
        'logg': [p[1] for p in params],
        'feh':  [p[2] for p in params],
        'vmic': [p[3] for p in params],
    }
    g = tm.grid(spectra, separable=(['vmic'] if sep_vmic else None), **param_dict)
    g.p0 = {'teff': teff, 'logg': logg, 'feh': feh, 'vmic': vmic}
    if not sep_vmic:
        g.interp_weights = np.asarray(weights, dtype=float)
    return g


def _apply_window(settings, xmin, xmax, pad):
    """Fill ``spectrum_params.aa_blue/aa_red`` from a wavelength window (unless a ``lam_file``
    is configured) and return the cache-path label. Mutates ``settings`` in place."""
    label = 'full'
    if xmin is not None and xmax is not None:
        lo, hi = (xmin, xmax) if xmin <= xmax else (xmax, xmin)
        sp = settings.setdefault('spectrum_params', {})
        if not sp.get('lam_file'):
            lo_i, hi_i = int(np.floor(lo)), int(np.ceil(hi))
            if pad is None:
                pad = min(max(0.5, 0.05 * (hi_i - lo_i)), 100)
            sp.setdefault('aa_blue', float(lo_i - pad))
            sp.setdefault('aa_red', float(hi_i + pad))
            label = f'w{lo_i}_{hi_i}'
    return label


def _dep_bucket(top, settings):
    """The ``dep_<hash>`` cache dir for an NLTE solve, keyed only on the departure-relevant
    settings so departures are reused across windows/line lists. ``None`` for LTE settings
    without the relevant sections."""
    dep_subset = {k: settings[k] for k in
                  ('atom_params', 'composition_params', 'opac_table_params', 'atmos_params',
                   'm3d_params', 'nlte_params', 'rt_iteration_params', 'rt_output_params',
                   'run_params')
                  if k in settings}
    dep_digest = hashlib.sha1(
        json.dumps(dep_subset, sort_keys=True, default=str).encode()).hexdigest()[:8]
    return Path(top) / f'dep_{dep_digest}'


def _interp_crop(grid, xmin, xmax, method):
    """Interpolate a node grid at its target point and crop to the window — the shared tail
    of ``compute_spectrum``/``compute_spectra``."""
    # a run can finish yet emit no spectrum (e.g. no R_sampling/daa/nnu/lam_file set, or
    # angle_quad='none'); flag it clearly instead of failing later on a missing attribute.
    no_spec = [s for s in grid.spectra if not getattr(s, 'has_spectrum', False)]
    if no_spec:
        folders = ', '.join(str(getattr(s, 'folder', '?')) for s in no_spec)
        raise RuntimeError(
            f"Synthesis produced no spectrum for: {folders}. Check the spectrum "
            f"settings (R_sampling / wavelength range) and the synthesis line list.")
    # pass LTE explicitly so an LTE run yields its (only) LTE spectrum without the "Cannot
    # get NLTE spectrum" warning Run.get_spectrum prints on the NLTE default
    lte = all(getattr(s, 'mode', 'LTE') == 'LTE' for s in grid.spectra)
    spec = grid.crop(xmin=xmin, xmax=xmax, method=method, LTE=lte)
    if len(spec) == 0:
        raise RuntimeError(
            f"Interpolated spectrum is empty over [{xmin}, {xmax}] Å — the synthesis "
            f"did not cover the requested window.")
    return spec


def _parse_abundances(value):
    """Parse the dispatch ``El=value`` abundance override into a stable mapping."""
    out = {}
    for term in str(value or '').replace(',', ' ').split():
        element, sep, raw = term.partition('=')
        if not sep or not element:
            raise ValueError(f"invalid abundance override {term!r}")
        key = element.lower().capitalize()
        if key in out:
            raise ValueError(f"duplicate abundance override for {key}")
        out[key] = float(raw)
    return out


def _format_abundances(values):
    return ' '.join(f'{element}={value:.12g}' for element, value in values.items())


PROVENANCE_NAME = 'provenance.json'


def _save_provenance(runs, nodes, settings):
    """Write what each finished run was asked for, beside its output.

    A run folder's name is rounded (Teff 5750 -> 't58'), so it cannot be read back as the
    point that was synthesised. The launcher holds the real values here -- they are what
    built the namelist -- so it writes them out rather than leaving the cache listing to
    recover them from formatted Fortran. Written for cache hits too, which back-fills a run
    that predates this file.
    """
    comp = (settings or {}).get('composition_params', {}) or {}
    try:
        items = list(_parse_abundances(comp.get('abundances', '')).items())
    except (TypeError, ValueError):
        items = []
    abund = {'items': items, 'relative': bool(comp.get('relative_abundances', False))}
    for jkey, teff_i, logg_i, feh_i, v, *_ in nodes:
        run = runs.get(jkey)
        if run is None:
            continue
        record = {'star': {'teff': float(teff_i), 'logg': float(logg_i),
                           'feh': float(feh_i), 'vmic': None if v is None else float(v)},
                  'abund': abund}
        try:
            (Path(run.folder) / PROVENANCE_NAME).write_text(json.dumps(record, indent=1))
        except OSError:
            # provenance is a convenience for the cache listing; a read-only or
            # vanished run directory must not fail the synthesis that produced it
            logging.exception('could not write %s for %s', PROVENANCE_NAME, run.folder)


def _foldable_ladder(settings_list, launch_kwargs, nlte):
    """Return ``(base_settings, manifest)`` when variants differ only in abundances."""
    variants = [copy.deepcopy(s or {}) for s in settings_list]
    abundance_maps = []
    comparable = []
    for settings in variants:
        comp = settings.setdefault('composition_params', {})
        try:
            abundance_maps.append(_parse_abundances(comp.pop('abundances', '')))
        except (TypeError, ValueError):
            return None
        comparable.append(settings)
    if any(item != comparable[0] for item in comparable[1:]):
        return None
    keys = set(abundance_maps[0])
    if not keys or any(set(values) != keys for values in abundance_maps[1:]):
        return None

    # The sampled table combines species and cannot apply an element-specific scale.
    effective = prepare_input(settings_list[0], launch_kwargs.get('default'))
    if any(bool((effective.get(section, {}) or {}).get('sample_linelist', False))
           for section in ('rt_iteration_params', 'rt_output_params')):
        return None

    base = copy.deepcopy(settings_list[0] or {})
    base_comp = base.setdefault('composition_params', {})
    base_comp['abundances'] = _format_abundances(abundance_maps[0])
    kind = 'nlte' if nlte else 'lte'
    manifest = [
        {'kind': 'continuum', 'offsets': {}},
        {'kind': kind, 'offsets': {}},
    ]
    if nlte:
        manifest.append({'kind': 'lte', 'offsets': {}})
    for values in abundance_maps:
        manifest.append({
            'kind': kind,
            'offsets': {element: values[element] - abundance_maps[0][element]
                        for element in abundance_maps[0]},
        })
    return base, manifest


class _LadderRunView:
    """One node spectrum as :class:`SpectrumGrid` expects it: set ``index`` of a folded run,
    or the run's own spectrum when ``index`` is None. ``output`` is the run's resolved
    folder as :func:`job_output` names it."""
    def __init__(self, run, index=None, output=None):
        self.run = run
        self.index = index
        self.has_spectrum = getattr(run, 'has_spectrum', False)
        self.mode = getattr(run, 'mode', 'LTE')
        self.folder = run.folder
        self.output = output

    @instance_cache
    def crop(self, **kwargs):
        """``run.crop`` of this node's spectrum; a workflow crops the same node to the same
        window at every trial star and lattice round, so the result is kept."""
        return self.run.crop(index=self.index, **kwargs)


class NodeStore:
    """The node spectra one workflow has synthesised on an abundance lattice.

    ``grid`` is a :class:`SpectrumGrid` over ``teff, logg, feh, vmic`` of the atmosphere
    node, ``snapshot`` (a snapshot node's tag, numbered by the store; 0 for a 1D node, whose
    vmic it replaces with 0) and ``abundance``, the element's value as written into the
    run's namelist ([X/Fe] or A(X), whichever the workflow uses). A node already in the grid is read from
    memory instead of being launched or read from disk again.

    Within one workflow a grid point has exactly one run. The pruning descriptor is
    therefore fixed by the first batch (``prune``), and a request whose run folder differs
    from the stored one for the same grid point raises instead of being served either.
    """
    def __init__(self):
        # a lookup table: no separable lattice is built, and it is never interpolated
        self.grid = tm.SpectrumGrid([], separable=[])
        self.prune = None
        self._pruned = False
        self._snapshots = {}
        self._selections = {}

    def select(self, star, atmos_files, **kwargs):
        """:func:`select_nodes` of ``star``, kept for the workflow, whose atmosphere grid
        is fixed: trial stars, probes and lattice rounds select the same stars repeatedly."""
        frozen = {k: tuple(sorted(v.items())) if isinstance(v, dict) else v
                  for k, v in kwargs.items()}
        key = (tuple(float(v) for v in star), tuple(sorted(frozen.items())))
        if key not in self._selections:
            self._selections[key] = select_nodes(star, atmos_files, **kwargs)
        return self._selections[key]

    def pruning(self, prune):
        """The workflow's pruning descriptor: ``prune`` on first use, the stored one after."""
        if not self._pruned:
            self.prune, self._pruned = prune, True
        return self.prune

    def point(self, node, abundance):
        jkey, teff, logg, feh, vmic = node[:5]
        snapshot = 0
        if vmic is None:
            snapshot = self._snapshots.setdefault(jkey[1], len(self._snapshots) + 1)
        return {'teff': teff, 'logg': logg, 'feh': feh,
                'vmic': 0.0 if vmic is None else vmic, 'snapshot': snapshot,
                'abundance': abundance}

    def get(self, node, abundance, output):
        """The stored view of ``node`` at ``abundance``, or None; raises when the stored
        run is not ``output``."""
        view = self.grid.find(**self.point(node, abundance))
        if view is not None and view.output != output:
            raise RuntimeError(
                f"Two runs for one grid point {self.point(node, abundance)}: "
                f"{view.output} and {output}.")
        return view

    def add(self, entries):
        """Store ``(node, abundance, view)`` entries in one grid update."""
        if not entries:
            return
        points = [self.point(node, abundance) for node, abundance, _view in entries]
        self.grid.add_data([view for _n, _a, view in entries],
                           **{k: [p[k] for p in points] for k in points[0]})


def _given_departures(nodes, dep_dirs):
    """The departure folders of ``nodes`` from a solved grid; every node must have one."""
    missing = [node[0] for node in nodes if node[0] not in dep_dirs]
    if missing:
        raise ValueError(f"The departure grid has no departures for nodes {missing}; solve "
                         "departures at this star.")
    return {node[0]: str(dep_dirs[node[0]]) for node in nodes}


def build_stellar_grid(top, star, atmos_files,
                       settings=None, absmet_files=None,
                       method='Delaunay', max_threads=None, launch_kwargs=None,
                       vmic_step=0.5, vmic_points=2, departures_top=None, nlte=None,
                       cancel_token=None, dep_dirs=None):
    """
    Build a tm.grid for the target stellar point by finding the enclosing
    model atmospheres, launching MULTI3D runs for them in parallel, and
    loading the results into a SpectrumGrid indexed by stellar parameters.

    vmic is *decoupled* from atmosphere selection: it barely affects the atmospheric
    structure but strongly affects line strengths, and MULTI3D applies it as a uniform
    namelist override. So atmospheres are bracketed only in (teff, logg, feh); for each
    selected node we take the single atmosphere nearest the query vmic and synthesise it
    at the two ``vmic_step``-spaced lattice points bracketing the query, then interpolate
    the spectra in vmic. This gives fine vmic sampling from a coarse atmosphere grid and
    lets vmic exceed the grid's range. (The selection itself lives in
    :func:`tumult.synthesis.atmos_select.select_nodes`.)

    Parameters
    ----------
    top : str or Path
        Output directory for the runs.
    star : tuple
        Target stellar parameters as (teff, logg, feh) or (teff, logg, feh, vmic).
    atmos_files : dict
        {(teff, logg, feh) or (teff, logg, feh, vmic): path} mapping. When keys carry a
        vmic, it is used only to pick the representative atmosphere per (teff, logg, feh).
    settings : dict, optional
        NameList settings (atom_params, nlte_params, …) merged into each
        per-atmosphere nml. atmos_params and composition_params are filled
        automatically per atmosphere.
    absmet_files : dict, optional
        {feh: path} mapping for absmet files.
    method : {'Delaunay', 'NaturalNeighbor'}
        How to select the enclosing atmosphere nodes (in teff/logg/feh).
    max_threads, launch_kwargs : optional
        Forwarded to _run_batch / tm.launch.
    vmic_step : float
        Microturbulence lattice spacing in km/s (default 0.5).
    vmic_points : int
        Number of vmic lattice points to synthesise per spatial node (default 2 = the
        linear bracket, historical behaviour). More points (3-4) cost proportionally
        more syntheses but let the returned grid interpolate vmic with monotone-cubic
        pchip instead of linearly: the vmic axis is flagged separable on the grid and
        ``interp_weights`` is not set (pchip weights depend on the data) — interpolate
        via ``grid.crop``, which ``compute_spectrum`` does automatically.
    departures_top : str or Path, optional
        When given *and* ``settings`` describe an NLTE solve, run NLTE in two steps: first
        solve departure coefficients (no line list, no output spectrum) for each atmosphere
        under ``departures_top``, saving them via the 'save' atmosphere format; then read
        them back in the line-synthesis step. ``departures_top`` is kept independent of the
        window/line list by the caller so departures are reused across them. When omitted,
        NLTE runs as a single combined solve+synthesis (legacy behaviour).
    dep_dirs : dict, optional
        ``{node key: departure folder}`` of departures solved earlier (a
        :class:`~tumult.departures.DepartureGrid`). The spectrum runs read them back and no
        departures are solved; every selected node must be present.
    nlte : bool, optional
        Whether this is an NLTE solve. Front-ends that carry the LTE/NLTE choice out of band
        (e.g. in namelist presets rather than ``settings``) must pass it explicitly, since it
        can't be recovered from ``settings`` alone. When ``None`` it is inferred from
        ``settings`` via :func:`_is_nlte` (works only when ``settings`` carries the
        iteration count, i.e. rt_iteration_params or legacy m3d_params).

    Returns
    -------
    tm.grid with a .target attribute {'teff', 'logg', 'feh', 'vmic'}.
    """
    launch_kwargs = launch_kwargs or {}
    if max_threads is None:
        max_threads = os.cpu_count() or 1

    nodes, vmic = select_nodes(star, atmos_files, method=method, vmic_step=vmic_step,
                               absmet_files=absmet_files, vmic_points=vmic_points)

    base = copy.deepcopy(settings or {})
    # explicit `nlte` wins (front-ends signal NLTE via presets, invisible to _is_nlte);
    # fall back to sniffing `settings` for direct library callers. Departures need a cache dir.
    if nlte is None:
        nlte = _is_nlte(base)
    nlte = departures_top is not None and bool(nlte)

    # NLTE step 1: solve departure coefficients per atmosphere (no line list, no output
    # spectrum), saved under `departures_top` (keyed independently of window/line list, so
    # reused across them). The synthesis step below reads them back.
    if dep_dirs is not None:
        dep_dirs, nlte = _given_departures(nodes, dep_dirs), True
    elif nlte:
        dep_runs = _run_batch(_dep_jobs(base, nodes, departures_top),
                              launch_kwargs, max_threads, cancel_token=cancel_token)
        dep_failed = [n[0] for n in nodes if dep_runs.get(n[0]) is None]
        if dep_failed:
            raise RuntimeError(f"Departure runs failed for atmospheres: {dep_failed}")
        dep_dirs = {n[0]: str(dep_runs[n[0]].folder) for n in nodes}
        _save_provenance(dep_runs, nodes, base)
    else:
        dep_dirs = {}

    jobs, job_keys, job_weights, job_params = _spec_jobs(base, nodes, top, nlte, dep_dirs)
    runs = _run_batch(jobs, launch_kwargs, max_threads, cancel_token=cancel_token)

    failed = [k for k in job_keys if runs.get(k) is None]
    if failed:
        raise RuntimeError(f"Runs failed for atmospheres: {failed}")
    _save_provenance(runs, nodes, base)

    return _assemble_grid(runs, star, vmic, job_keys, job_weights, job_params, vmic_points)


def compute_spectrum(top, star, atmos_files, *, settings=None, absmet_files=None,
                     method='Delaunay', max_threads=None, launch_kwargs=None,
                     xmin=None, xmax=None, pad=None, vmic_step=0.5, vmic_points=2,
                     nlte=None, cancel_token=None, return_runs=False, dep_dirs=None):
    """Compute a single spectrum at an arbitrary stellar point (Teff, logg, FeH, vmic).

    Thin wrapper over :func:`build_stellar_grid`: it synthesises the enclosing model
    atmospheres in parallel and interpolates the result at ``star``. Whether the run is NLTE
    is taken from the explicit ``nlte`` flag when given (front-ends that carry LTE/NLTE in
    namelist presets must pass it — it can't be recovered from ``settings``), else inferred
    from ``settings`` via :func:`_is_nlte`. NLTE runs use the two-step departures path.
    This helper also fills in the synthesis wavelength window from ``xmin``/``xmax`` so the
    Fortran code computes just the visible range.

    Parameters
    ----------
    top, star, atmos_files, settings, absmet_files, method, max_threads, launch_kwargs
        Forwarded to :func:`build_stellar_grid` (``settings`` is the fully assembled nml
        dict, keyed by namelist section).
    xmin, xmax : float, optional
        Wavelength window in Angstrom. Sets ``spectrum_params.aa_blue/aa_red`` (unless a
        ``lam_file`` is configured) and crops the interpolated spectrum to this range.
    pad : float, optional
        Extra Angstrom synthesised on each side of [xmin, xmax] so the cropped window has
        clean edges for convolution. Defaults to ``max(0.5, 0.05 * (xmax - xmin))``.
    vmic_points : int, optional
        vmic lattice points per spatial node (default 2 = linear bracket). With 3-4 the
        vmic axis is interpolated by monotone-cubic pchip via the grid's separable-axis
        machinery — more accurate for saturated lines, at the cost of proportionally
        more syntheses (mitigated by the run cache: lattice levels are shared across
        queries). In that mode the spatial weights are re-derived on the node stencil
        by ``grid.crop`` (identical for Delaunay; NaturalNeighbor weights are then
        computed on the stencil instead of the full atmosphere grid).

    return_runs : bool
        Also return the producing run directories, for consuming run artifacts.
    dep_dirs : dict, optional
        Departures solved earlier, per node (see :func:`build_stellar_grid`).

    Returns
    -------
    Spectrum or (Spectrum, tuple[Path, ...])
        The continuum-normalised spectrum, optionally with its producing runs.
    """
    settings = copy.deepcopy(settings or {})
    top = Path(top)

    label = _apply_window(settings, xmin, xmax, pad)

    # build_stellar_grid names cached runs by stellar params only, so a run computed with
    # different settings (window, R, abundances, mode, line list, ...) would be silently
    # reused. Key the output path on a hash of the settings so any change re-synthesises;
    # identical settings still hit the disk cache.
    digest = hashlib.sha1(
        json.dumps(settings, sort_keys=True, default=str).encode()).hexdigest()[:8]
    syn_top = top / f'{label}_{digest}'

    # explicit `nlte` wins (front-ends signal NLTE via presets, invisible to _is_nlte);
    # fall back to sniffing `settings` for direct library callers.
    if nlte is None:
        nlte = _is_nlte(settings)

    # NLTE departure runs are independent of the window and line list, so key their cache on
    # only the departure-relevant settings and place them in a sibling `dep_<hash>` bucket —
    # reused across windows/line lists (see build_stellar_grid's two-step path).
    departures_top = _dep_bucket(top, settings) if nlte else None

    grid = build_stellar_grid(str(syn_top), star, atmos_files, settings=settings,
                              absmet_files=absmet_files, method=method,
                              max_threads=max_threads, launch_kwargs=launch_kwargs,
                              vmic_step=vmic_step, vmic_points=vmic_points,
                              departures_top=str(departures_top) if departures_top else None,
                              nlte=nlte, cancel_token=cancel_token, dep_dirs=dep_dirs)

    spec = _interp_crop(grid, xmin, xmax, method)
    if return_runs:
        return spec, tuple(Path(run.folder).resolve() for run in grid.spectra)
    return spec


def normalise_windows(windows, pad=None):
    """Validate, pad, sort and merge synthesis windows.

    ``windows`` is a sequence of ``(blue, red)`` bounds in Angstrom. Each is padded so a
    cropped window has clean edges for convolution, then overlapping or touching windows
    are merged: they would otherwise produce segments that share grid points, and the
    Fortran side merges them anyway. Returns the merged ``[(blue, red), ...]`` sorted by
    blue edge.

    The padding rule matches the single-window path (:func:`_apply_window`), so one window
    through here and the same window through ``compute_spectrum`` cover the same range.
    """
    if not windows:
        raise ValueError("compute_windows needs at least one (blue, red) window")
    padded = []
    for w in windows:
        blue, red = (float(w[0]), float(w[1]))
        if not (np.isfinite(blue) and np.isfinite(red)):
            raise ValueError(f"window {w} has a non-finite bound")
        if blue > red:
            blue, red = red, blue
        if blue <= 0:
            raise ValueError(f"window {w} has a non-positive wavelength")
        lo_i, hi_i = int(np.floor(blue)), int(np.ceil(red))
        p = pad
        if p is None:
            p = min(max(0.5, 0.05 * (hi_i - lo_i)), 100)
        padded.append((float(lo_i - p), float(hi_i + p)))
    padded.sort()
    merged = [padded[0]]
    for blue, red in padded[1:]:
        if blue <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], red))
        else:
            merged.append((blue, red))
    return merged


def window_file_text(windows):
    """The ``window_file`` body for ``windows``, as ``spectrum_params`` reads it.

    Deterministic, so its content can key the run cache: two requests for the same
    windows must reuse the same synthesis, and a changed window must not.
    """
    lines = ["# blue_AA  red_AA"]
    lines += [f"{blue:.4f}  {red:.4f}" for blue, red in windows]
    return "\n".join(lines) + "\n"


def write_window_file(folder, windows):
    """Write the ``window_file`` of normalised ``windows`` under ``folder``; its path.

    The file is named by its content, so the path a namelist records identifies the
    windows and two runs over the same windows share one cache entry.
    """
    text = window_file_text(windows)
    path = Path(folder) / f"{hashlib.sha1(text.encode()).hexdigest()[:8]}.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path.resolve()


def compute_windows(top, star, atmos_files, *, windows, R, settings=None, absmet_files=None,
                    method='Delaunay', max_threads=None, launch_kwargs=None, pad=None,
                    vmic_step=0.5, vmic_points=2, nlte=None, cancel_token=None):
    """Synthesise several disjoint wavelength windows in one run, one spectrum each.

    The windows are normalised (:func:`normalise_windows`) and written to a window file
    the run reads, so the Fortran side builds one constant-``R`` segment per window and
    splits the output patches at window boundaries. Overlapping requests merge into one segment; the
    result is still cropped back to every window that was asked for, so the caller gets a
    list aligned with ``windows`` whether or not any merging happened.

    Returns ordinary contiguous spectra, one per requested window, in the order given.
    Each is what synthesising that window on its own would have produced.

    ``R`` is required: a window file fixes the sampling per window, and mixing that with
    ``daa``/``nnu`` would leave the grid ambiguous.
    """
    if R is None or R <= 0:
        raise ValueError("compute_windows needs a positive R (constant-resolution grid)")
    top = Path(top)
    settings = copy.deepcopy(settings or {})
    merged = normalise_windows(windows, pad)

    # The window file is a generated input, so the cache must key on its CONTENT: a path
    # alone would let two different window sets share a digest.
    text = window_file_text(merged)
    sp = settings.setdefault('spectrum_params', {})
    sp['R_sampling'] = float(R)
    digest = hashlib.sha1(json.dumps(
        [settings, text], sort_keys=True, default=str).encode()).hexdigest()[:8]
    window_path = (Path((launch_kwargs or {}).get('cwd') or '.') / top /
                   '_windows' / f'{digest}.txt').resolve()
    window_path.parent.mkdir(parents=True, exist_ok=True)
    window_path.write_text(text)
    sp['window_file'] = str(window_path)

    if nlte is None:
        nlte = _is_nlte(settings)
    syn_top = top / f'win{len(merged)}_{digest}'
    departures_top = _dep_bucket(top, settings) if nlte else None

    grid = build_stellar_grid(str(syn_top), star, atmos_files, settings=settings,
                              absmet_files=absmet_files, method=method,
                              max_threads=max_threads, launch_kwargs=launch_kwargs,
                              vmic_step=vmic_step, vmic_points=vmic_points,
                              departures_top=str(departures_top) if departures_top else None,
                              nlte=nlte, cancel_token=cancel_token)

    # One spectrum per REQUESTED window, cropped out of the segment covering it. Cropping
    # inside a segment never crosses a gap, so these are ordinary contiguous spectra.
    return [_interp_crop(grid, float(min(w)), float(max(w)), method) for w in windows]


def _select_all(stars, atmos_files, store=None, **kwargs):
    """Node selections for every star in ``stars`` and the union of their nodes, in order."""
    select = select_nodes if store is None else store.select
    selections = [select(st, atmos_files, **kwargs) for st in stars]
    union = {}
    for nodes, _vmic in selections:
        for node in nodes:
            union.setdefault(node[0], node)
    return selections, list(union.values())


def _run_plans(plans, nodes, selections, stars, *, top, launch_kwargs, max_threads,
               cancel_token, given, store, abundances, set_index, crops, method,
               vmic_points):
    """Run every (plan x node) pair a store does not hold yet and interpolate each star.

    ``plans`` are ``(members, base, syn_top, nlte, node_lists)``: one run per node per plan,
    its ``members`` the variants it carries (several only for a folded ladder, whose
    variant ``slot`` is read from set ``set_index(nlte, slot)``; None reads a run's own
    spectrum). ``nodes`` is the union of the node ``selections`` of ``stars``. Departures
    come from ``given`` (per node) or are solved in one batch before the spectra.

    Returns ``(results, views, dep_runs)`` with ``results[star][variant][crop]`` for the
    ``(xmin, xmax)`` pairs in ``crops``.
    """
    stored = {}
    if store is not None:
        for gi, (members, base, syn_top, _nlte, _lists) in enumerate(plans):
            abundance = abundances[members[0]]
            for node in nodes:
                name = syn_top / _node_name(node[0], *node[1:5])
                view = store.get(node, abundance, job_output(name, base, launch_kwargs))
                if view is not None:
                    stored[(gi, node[0])] = view
    missing = [[node for node in nodes if (gi, node[0]) not in stored]
               for gi in range(len(plans))]

    # every plan's departures in one batch, then every plan's spectra in one more, so a
    # ladder fans (rung x node) across the pool instead of a rung at a time
    dep_dirs, dep_runs = {}, {}
    if given is not None:
        dep_dirs = {(gi, key): folder for gi in range(len(plans))
                    for key, folder in given.items()}
    dep_jobs = [((gi, jkey), path, nml)
                for gi, (_m, base, _top, nlte, _lists) in enumerate(plans)
                if nlte and given is None
                for jkey, path, nml in _dep_jobs(base, missing[gi], _dep_bucket(top, base))]
    if dep_jobs:
        dep_runs = _run_batch(dep_jobs, launch_kwargs, max_threads,
                              cancel_token=cancel_token)
        failed = [key for key, _path, _nml in dep_jobs if dep_runs.get(key) is None]
        if failed:
            raise RuntimeError(f"Departure runs failed: {failed}")
        dep_dirs = {key: str(dep_runs[key].folder) for key, _path, _nml in dep_jobs}

    spec_jobs = []
    for gi, (_members, base, syn_top, nlte, node_lists) in enumerate(plans):
        dep_v = {n[0]: dep_dirs[(gi, n[0])] for n in missing[gi]} if nlte else {}
        jobs = _spec_jobs(base, missing[gi], syn_top, nlte, dep_v, node_lists)[0]
        spec_jobs += [((gi, jkey), path, nml) for jkey, path, nml in jobs]
    runs = _run_batch(spec_jobs, launch_kwargs, max_threads, cancel_token=cancel_token)
    failed = [key for key, _path, _nml in spec_jobs if runs.get(key) is None]
    if failed:
        raise RuntimeError(f"Runs failed: {failed}")

    n_variants = sum(len(plan[0]) for plan in plans)
    views = {}           # (variant, node key) -> view
    new_entries = []
    for gi, (members, _base, _top, nlte, _lists) in enumerate(plans):
        outputs = {key: job_output(path, nml, launch_kwargs)
                   for (g, key), path, nml in spec_jobs if g == gi}
        for slot, vi in enumerate(members):
            for node in nodes:
                key = node[0]
                view = stored.get((gi, key))
                if view is None:
                    view = _LadderRunView(runs[(gi, key)], set_index(nlte, slot),
                                          outputs[key])
                    if store is not None:
                        new_entries.append((node, abundances[vi], view))
                views[(vi, key)] = view
    if store is not None:
        store.add(new_entries)

    results = []
    for (sel_nodes, vmic), st in zip(selections, stars):
        keys = [node[0] for node in sel_nodes]
        weights = [node[5] for node in sel_nodes]
        params = [(node[1], node[2], node[3], 0.0 if node[4] is None else node[4])
                  for node in sel_nodes]
        per_variant = []
        for vi in range(n_variants):
            grid = _assemble_grid({key: views[(vi, key)] for key in keys}, st, vmic, keys,
                                  weights, params, vmic_points)
            per_variant.append([_interp_crop(grid, lo, hi, method) for lo, hi in crops])
        results.append(per_variant)
    return results, list(views.values()), dep_runs


def compute_window_batch(top, star, atmos_files, *, settings_list, windows, R,
                         absmet_files=None, method='Delaunay', max_threads=None,
                         launch_kwargs=None, pad=None, vmic_step=0.5, vmic_points=2,
                         nlte=None, cancel_token=None, prune=None, fold=False,
                         preparation=None, return_runs=False, dep_dirs=None,
                         store=None, abundances=None, probes=()):
    """Both axes in one run: an abundance ladder across several wavelength windows.

    ``settings_list`` is a ladder in the sense of :func:`compute_spectra` — variants
    differing only in one element's abundance — and ``windows`` a sequence of
    ``(blue, red)`` bounds. The ladder folds onto the run's set axis and the windows onto
    its output patches, which are orthogonal, so a single synthesis covers every
    (variant, window) pair: the atmosphere, the equation of state, the background and, in
    NLTE, the statistical-equilibrium solve are each paid once instead of once per pair.

    Returns ``results[variant][window]``, each an ordinary contiguous spectrum matching
    the requested window order.

    Folding solves NLTE departures at the first variant's abundance and scales the others,
    which is not trusted for NLTE and is therefore an explicit opt-in (``fold=True``).

    Falls back to :func:`compute_spectra` per window when the variants are not a foldable
    ladder (different line lists, elements switched on and off, ``sample_linelist``), so a
    caller does not have to know which of its requests can fold.

    ``prune`` (optional) asks for the per-node line-list pruning described in
    :func:`_node_linelists`, and is part of these runs' cache identity.

    ``dep_dirs`` supplies departures solved earlier, per node (see
    :func:`build_stellar_grid`); no departures are then solved.

    ``fold=False`` (the default) keeps the windows on the patches but gives every variant its own run,
    scheduled together with the others: one :func:`_run_batch` over every
    ``(variant × node)`` pair instead of one per variant. That trades the folded run's
    shared atmosphere and background for a cache key that is one abundance rather than a
    whole ladder, which is what lets a rung on an abundance lattice be reused by a later
    solver iteration that does not repeat the ladder around it. Each variant's runs are
    byte-identical to the ones it would get from its own ``fold=True`` call with that one
    variant, so the two schedulings share a cache.

    ``store`` (a :class:`NodeStore`, with ``fold=False``) holds the node runs of earlier
    calls in the same workflow; ``abundances`` gives each variant's lattice coordinate in
    it. Only the (variant x node) runs missing from the store are launched, together with
    their departure solves, and the new ones are added to it. The store also fixes the
    pruning descriptor for the whole workflow, overriding ``prune``.

    ``probes`` are further stars interpolated from the same variants, their nodes run
    together with the target's; the results then carry one entry per star,
    ``results[star][variant][window]``, the target first.
    """
    if R is None or R <= 0:
        raise ValueError("compute_window_batch needs a positive R")
    if not settings_list:
        return []
    top = Path(top)
    launch_kwargs = launch_kwargs or {}
    if max_threads is None:
        max_threads = os.cpu_count() or 1
    merged = normalise_windows(windows, pad)
    if store is not None:
        if fold:
            raise ValueError("A node store holds one run per abundance; use fold=False.")
        if len(abundances or ()) != len(settings_list):
            raise ValueError("A node store needs one abundance per variant.")
        prune = store.pruning(prune)

    stars = [star, *probes]
    selections, nodes = _select_all(stars, atmos_files, store, method=method,
                                    vmic_step=vmic_step,
                                    absmet_files=absmet_files, vmic_points=vmic_points)

    text = window_file_text(merged)
    prepared = []
    for settings in settings_list:
        item = copy.deepcopy(settings or {})
        sp = item.setdefault('spectrum_params', {})
        sp['R_sampling'] = float(R)
        prepared.append(item)
    nlte_values = [bool(_is_nlte(s) if nlte is None else nlte) for s in prepared]

    # one group of every variant when folding, else one group per variant: the rest of the
    # path is the same, and a group of one is exactly what a single-variant call produces
    groups = [list(range(len(prepared)))] if fold else [[i] for i in range(len(prepared))]
    folds = []
    for members in groups:
        nlte_group = {nlte_values[i] for i in members}
        folds.append(_foldable_ladder([prepared[i] for i in members], launch_kwargs,
                                      nlte_values[members[0]])
                     if len(nlte_group) == 1 else None)
    if dep_dirs is not None:
        dep_dirs = _given_departures(nodes, dep_dirs)
        if not all(nlte_values):
            raise ValueError("Solved grid departures are read back by NLTE settings only.")
    if any(item is None for item in folds):
        if store is not None or probes:
            raise ValueError("A node store and probe stars need variants that differ only "
                             "in abundances.")
        # not a ladder: nothing to fold onto the set axis, so fall back to the existing
        # per-window batching rather than pretending this is one run
        per_window, folders = [], []
        for w in windows:
            specs, runs_w = compute_spectra(
                top, star, atmos_files, settings_list=settings_list,
                absmet_files=absmet_files, method=method, max_threads=max_threads,
                launch_kwargs=launch_kwargs, xmin=float(min(w)), xmax=float(max(w)),
                pad=pad, vmic_step=vmic_step, vmic_points=vmic_points, nlte=nlte,
                cancel_token=cancel_token, prune=prune, preparation=preparation,
                return_runs=True, dep_dirs=dep_dirs)
            per_window.append(specs)
            folders.extend(runs_w)
        values = [list(specs) for specs in zip(*per_window)]
        return (values, tuple(dict.fromkeys(folders))) if return_runs else values

    plans = []      # per group: (members, base, syn_top, nlte_v, node_lists)
    for members, (base, manifest) in zip(groups, folds):
        manifest_text = json.dumps(manifest, sort_keys=True, separators=(',', ':'))
        node_lists = _node_linelists(top, base, nodes, merged, prune, launch_kwargs,
                                     preparation, cancel_token)
        # the pruning only enters the digest when there is one, so an unpruned batch keeps the
        # cache identity (and the folders) it had before per-node pruning existed
        keys = [base, manifest_text, text]
        if node_lists:      # a prune or a prepared list: the selection enters the digest
            keys.append(str(preparation.folder) if preparation is not None else
                        _prune_digest(base, merged, prune))
        digest = hashlib.sha1(
            json.dumps(keys, sort_keys=True, default=str).encode()).hexdigest()[:8]
        manifest_path = (Path((launch_kwargs or {}).get('cwd') or '.') / top /
                         '_set_manifests' / f'{digest}.txt').resolve()
        write_set_manifest(manifest_path, manifest)
        window_path = (Path((launch_kwargs or {}).get('cwd') or '.') / top /
                       '_windows' / f'{digest}.txt').resolve()
        window_path.parent.mkdir(parents=True, exist_ok=True)
        window_path.write_text(text)
        base.setdefault('output_params', {})['set_file'] = str(manifest_path)
        base.setdefault('spectrum_params', {})['window_file'] = str(window_path)
        plans.append((members, base, top / f'win{len(merged)}_{digest}',
                      nlte_values[members[0]], node_lists))

    results, views, _dep_runs = _run_plans(
        plans, nodes, selections, stars, top=top, launch_kwargs=launch_kwargs,
        max_threads=max_threads, cancel_token=cancel_token, given=dep_dirs, store=store,
        abundances=abundances, set_index=lambda nlte, slot: (3 if nlte else 2) + slot,
        crops=[(float(min(w)), float(max(w))) for w in windows], method=method,
        vmic_points=vmic_points)
    results = results if probes else results[0]
    if return_runs:
        folders = tuple(dict.fromkeys(str(Path(view.folder).resolve()) for view in views))
        return results, folders
    return results


def compute_spectra(top, star, atmos_files, *, settings_list, absmet_files=None,
                    method='Delaunay', max_threads=None, launch_kwargs=None,
                    xmin=None, xmax=None, pad=None, vmic_step=0.5, vmic_points=2,
                    nlte=None, cancel_token=None, prune=None, preparation=None,
                    return_runs=False, fold=False, dep_dirs=None, store=None,
                    abundances=None):
    """Batched :func:`compute_spectrum` for a list of settings variants that share the same
    star/atmospheres (e.g. a curve-of-growth's abundance ladder).

    Every ``(atmosphere-node × variant)`` dispatch run uses the shared execution
    service: one :func:`_run_batch` for the NLTE departure solves and one for the
    spectrum solves (just one for LTE). Atmosphere/vmic node selection is done once and shared
    across variants. Returns ``[Spectrum]`` aligned with ``settings_list``.

    Each variant's runs are byte-identical to the ones :func:`compute_spectrum` would issue
    for that settings dict — same window, same hash-keyed cache folders, same namelists — so
    the batch is a pure scheduling change: results match the per-variant serial calls exactly
    and the disk cache is shared with them. Departures are solved per variant (the dep cache
    key includes ``composition_params``), so an abundance ladder gets a full NLTE solve at
    every point, only run concurrently instead of one after another.

    ``prune`` (optional) asks for the per-node line-list pruning described in
    :func:`_node_linelists`, and is part of these runs' cache identity.

    Variants differing only in abundances fold into one run whose departures are solved at
    the first variant's abundance and whose other rungs are opacity scalings of it (see
    :func:`_foldable_ladder`). That is an explicit opt-in (``fold=True``) for LTE callers:
    by default every variant keeps its own run, so an NLTE ladder gets a departure solve
    per rung.

    ``dep_dirs`` supplies departures solved earlier, per node (see
    :func:`build_stellar_grid`); no departures are then solved.

    ``store`` and ``abundances`` keep the unfolded variants' runs for a workflow, as in
    :func:`compute_window_batch`.
    """
    if not settings_list:
        return []
    top = Path(top)
    launch_kwargs = launch_kwargs or {}
    if max_threads is None:
        max_threads = os.cpu_count() or 1
    if store is not None:
        if fold:
            raise ValueError("A node store holds one run per abundance; use fold=False.")
        if len(abundances or ()) != len(settings_list):
            raise ValueError("A node store needs one abundance per variant.")
        prune = store.pruning(prune)

    # atmosphere/vmic nodes depend only on the star + grid, not the per-variant settings —
    # select once and share across variants.
    nodes, vmic = select_nodes(star, atmos_files, method=method, vmic_step=vmic_step,
                               absmet_files=absmet_files, vmic_points=vmic_points)

    # Apply the window before comparing variants: it is part of both the fold condition
    # and cache identity.
    windowed = []
    labels = []
    for settings in settings_list:
        item = copy.deepcopy(settings or {})
        labels.append(_apply_window(item, xmin, xmax, pad))
        windowed.append(item)
    nlte_values = [bool(_is_nlte(s) if nlte is None else nlte) for s in windowed]

    given = None if dep_dirs is None else _given_departures(nodes, dep_dirs)
    if given is not None:
        nlte_values = [True] * len(windowed)
        fold = False
    fold = (_foldable_ladder(windowed, launch_kwargs, nlte_values[0])
            if fold and len(set(nlte_values)) == 1 and len(set(labels)) == 1 else None)
    prune_windows = ([(float(min(xmin, xmax)), float(max(xmin, xmax)))]
                     if xmin is not None and xmax is not None else [])
    if fold is not None:
        base, manifest = fold
        manifest_text = json.dumps(manifest, sort_keys=True, separators=(',', ':'))
        node_lists = _node_linelists(top, base, nodes, prune_windows, prune, launch_kwargs,
                                     preparation, cancel_token) \
            if prune_windows else None
        keys = [base, manifest_text]
        if node_lists:      # a prune or a prepared list: the selection enters the digest
            keys.append(str(preparation.folder) if preparation is not None else
                        _prune_digest(base, prune_windows, prune))
        digest = hashlib.sha1(
            json.dumps(keys, sort_keys=True, default=str).encode()).hexdigest()[:8]
        manifest_path = (Path((launch_kwargs or {}).get('cwd') or '.') / top /
                     '_set_manifests' / f'{digest}.txt').resolve()
        write_set_manifest(manifest_path, manifest)
        base.setdefault('output_params', {})['set_file'] = str(manifest_path)
        syn_top = top / f'{labels[0]}_{digest}'
        nlte_v = nlte_values[0]
        dep_top = _dep_bucket(top, base) if nlte_v else None

        dep_dirs = {}
        if nlte_v:
            dep_jobs = _dep_jobs(base, nodes, dep_top)
            dep_runs = _run_batch(dep_jobs, launch_kwargs, max_threads,
                                  cancel_token=cancel_token)
            failed = [key for key, _path, _nml in dep_jobs if dep_runs.get(key) is None]
            if failed:
                raise RuntimeError(f"Departure runs failed: {failed}")
            dep_dirs = {key: str(dep_runs[key].folder) for key, _path, _nml in dep_jobs}

        jobs, job_keys, job_weights, job_params = _spec_jobs(
            base, nodes, syn_top, nlte_v, dep_dirs, node_lists)
        runs = _run_batch(jobs, launch_kwargs, max_threads, cancel_token=cancel_token)
        failed = [key for key in job_keys if runs.get(key) is None]
        if failed:
            raise RuntimeError(f"Runs failed: {failed}")

        first_ladder = 3 if nlte_v else 2
        results = []
        for vi in range(len(settings_list)):
            views = {key: _LadderRunView(runs[key], first_ladder + vi)
                     for key in job_keys}
            grid = _assemble_grid(views, star, vmic, job_keys, job_weights,
                                  job_params, vmic_points)
            results.append(_interp_crop(grid, xmin, xmax, method))
        return results

    # Per variant: assemble window/hash-keyed tops exactly as compute_spectrum does.
    plans = []      # (members, settings, syn_top, nlte_v, node_linelists)
    for vi, (s, label, nlte_v) in enumerate(zip(windowed, labels, nlte_values)):
        node_lists = _node_linelists(top, s, nodes, prune_windows, prune, launch_kwargs,
                                     preparation, cancel_token) \
            if prune_windows else None
        prep_key = (str(preparation.folder) if preparation is not None else
                    _prune_digest(s, prune_windows, prune))
        keys = ([s, prep_key] if node_lists is not None
                else s)
        digest = hashlib.sha1(
            json.dumps(keys, sort_keys=True, default=str).encode()).hexdigest()[:8]
        plans.append(([vi], s, top / f'{label}_{digest}', bool(nlte_v), node_lists))

    results, views, dep_runs = _run_plans(
        plans, nodes, [(nodes, vmic)], [star], top=top, launch_kwargs=launch_kwargs,
        max_threads=max_threads, cancel_token=cancel_token, given=given, store=store,
        abundances=abundances, set_index=lambda nlte, slot: None, crops=[(xmin, xmax)],
        method=method, vmic_points=vmic_points)
    results = [specs[0] for specs in results[0]]
    if return_runs:
        folders = tuple(dict.fromkeys(
            str(Path(item.folder).resolve())
            for item in (*dep_runs.values(), *views)))
        return results, folders
    return results


def adaptive_fit(target, synth_point_fn, axes, guesses, *, post_fit=(), spacings=None,
                 axis_nodes=None, atmos_bounds=None, max_rounds=6, tol=None, xtol=None,
                 fit_kwargs=None, separable=None, cancel_token=None):
    """Fit ``target`` by adaptive grid synthesis.

    Builds a minimal grid (2 points per fitted axis, bracketing the guess), synthesises each
    point via ``synth_point_fn``, and fits every axis at once with :meth:`SpectrumGrid.fit`
    (which also fits the continuous ``post_fit`` params — vshift/vrot/zeta_rt — exactly, with
    no grid axis). It then expands any axis whose optimum sits on a boundary, or subdivides
    around an interior optimum that is still coarsely sampled, and refits — until every axis
    optimum is interior and resolved to ``tol`` (or ``max_rounds`` is reached).

    Parameters
    ----------
    target : Spectrum
        Observed spectrum to match (compared on the synthetic grid via ``crop(newx=...)``).
    synth_point_fn : callable
        ``synth_point_fn(point: dict) -> Spectrum`` synthesising a normalised spectrum at one
        grid point (keys are ``axes``). May return ``None`` for an un-synthesisable point.
        Memoised here on the rounded point, so repeated points cost nothing.
    axes : sequence of str
        Names of the fitted grid dimensions (e.g. 'teff', 'logg', 'feh', 'vmic', 'Mg').
    guesses : dict
        Initial value for every fitted name — grid ``axes`` *and* ``post_fit`` params.
    post_fit : sequence of str, optional
        Continuous post-processing params (vshift/vrot/zeta_rt) fit by ``Spectrum.fit``.
    spacings : dict, optional
        Step for each *continuous* axis (vmic, abundances) — used to bracket the guess and
        expand outward, with sub-step refinement around an interior optimum. Default 1.0.
    axis_nodes : dict, optional
        ``{axis: sorted values}`` for *discrete* axes whose sample points are fixed (the
        atmosphere-grid nodes for teff/logg/feh). These bracket the guess and expand to
        adjacent nodes but are never sub-sampled (finer than a node is just interpolation
        between the same two syntheses).
    atmos_bounds : dict, optional
        ``{axis: (lo, hi)}`` hard clamp for continuous-axis expansion (``None`` = unbounded),
        e.g. ``(0, None)`` for vmic. Discrete axes are bounded by their node list.
    max_rounds : int, optional
    tol : dict, optional
        Per-axis convergence tolerance for grid expansion/subdivision (default:
        ``spacings[axis] / 4``). This governs how finely the *grid* is sampled, not the
        optimiser precision (see ``xtol``).
    xtol : dict, optional
        Per-parameter absolute optimiser tolerance, forwarded to ``Spectrum.fit`` (e.g.
        ``{"Mg": 0.01, "vshift": 0.1}``). Also the precision the returned best-fit values are
        rounded to. Velocities default to 0.1; everything else to 0.01.
    fit_kwargs : dict, optional
        Extra kwargs forwarded to ``grid.fit`` (e.g. ``method``, ``integrate``, ``sigma``).
    separable : sequence of str, optional
        Axis names to interpolate with monotone-cubic pchip instead of linearly (the
        grid's separable-axis machinery; names are whatever you called the axes).
        Physically right for abundance-like axes, where the flux responds monotonically
        and the curve of growth is strongly curved. Requires every sampled grid point
        to synthesise; if one fails (``synth_point_fn`` returned ``None``) the round
        falls back to linear interpolation with a printed note.

    Returns
    -------
    (best : dict, fitted : Spectrum, grid : SpectrumGrid)
    """
    axes = list(axes)
    post_fit = list(post_fit)
    spacings = dict(spacings or {})
    axis_nodes = {k: sorted(float(x) for x in v) for k, v in (axis_nodes or {}).items()}
    atmos_bounds = dict(atmos_bounds or {})
    fit_kwargs = dict(fit_kwargs or {})
    if xtol:
        fit_kwargs.setdefault("xtol", dict(xtol))   # optimiser precision per parameter
    tol = dict(tol or {})

    def _clamp(axis, v):
        lo, hi = atmos_bounds.get(axis, (None, None))
        if lo is not None:
            v = max(v, lo)
        if hi is not None:
            v = min(v, hi)
        return v

    def _bracket_nodes(nodes, g):
        """The two nodes bracketing ``g`` (or the nearest two if ``g`` is off the ends)."""
        i = int(np.searchsorted(nodes, g))
        i = max(1, min(i, len(nodes) - 1))
        return [nodes[i - 1], nodes[i]]

    cache = {}

    def _synth(point):
        key = tuple(sorted((k, round(float(v), 6)) for k, v in point.items()))
        if key not in cache:
            cache[key] = synth_point_fn(dict(point))
        return cache[key]

    # initial samples per axis: discrete axes bracket the guess at their nodes; continuous
    # axes take 2 points bracketing the guess (re-spread if clamping collapses them).
    samples = {}
    for a in axes:
        g = float(guesses[a])
        if a in axis_nodes:
            nodes = axis_nodes[a]
            samples[a] = _bracket_nodes(nodes, g) if len(nodes) >= 2 else list(nodes) or [g]
        else:
            s = float(spacings.get(a, 1.0))
            pts = sorted({_clamp(a, g - s / 2), _clamp(a, g + s / 2)})
            if len(pts) < 2:
                base = pts[0] if pts else g
                pts = sorted({base, _clamp(a, base + s), _clamp(a, base - s)})
            samples[a] = pts

    best = fitted = grid = None
    for _round in range(max_rounds):
        if cancel_token is not None:
            cancel_token.check()   # abort promptly between refinement rounds
        # synthesise the Cartesian product of the current samples
        spectra, params = [], {a: [] for a in axes}
        for combo in itertools.product(*(samples[a] for a in axes)):
            point = dict(zip(axes, combo))
            spec = _synth(point)
            if spec is None:
                continue
            spectra.append(spec)
            for a, v in point.items():
                params[a].append(v)
        # need >=2 points to interpolate a grid axis; a post-params-only fit (no axes) is
        # a single synthesis that Spectrum.fit transforms in place.
        if len(spectra) < (2 if axes else 1):
            raise RuntimeError("adaptive_fit: too few synthesised grid points to fit.")

        if separable is not None:
            try:
                grid = tm.grid(spectra, separable=list(separable), **params)
            except ValueError as err:
                print(f"adaptive_fit: separable interpolation unavailable this round "
                      f"({err}); falling back to linear.")
                grid = tm.grid(spectra, **params)
        else:
            grid = tm.grid(spectra, **params)
        p0 = {k: (list(guesses[k]) if k == "continuum" else float(guesses[k]))
              for k in axes + post_fit if k in guesses}
        best, fitted = grid.fit(target, fit=axes + post_fit, p0=p0, **fit_kwargs)

        # per-axis: expand at a boundary optimum; for continuous axes also subdivide a coarse
        # interior bracket. Discrete (node) axes only expand to the next node.
        changed = False
        for a in axes:
            pts = samples[a]
            lo, hi = pts[0], pts[-1]
            v = float(best[a])
            if a in axis_nodes:
                nodes = axis_nodes[a]
                t = float(tol.get(a, 0.1 * (hi - lo if hi > lo else 1.0)))
                if v <= lo + t:                       # extend to the previous node
                    i = int(np.searchsorted(nodes, lo))
                    if i > 0 and nodes[i - 1] not in pts:
                        pts.insert(0, nodes[i - 1])
                        changed = True
                elif v >= hi - t:                     # extend to the next node
                    i = int(np.searchsorted(nodes, hi))
                    if i < len(nodes) - 1 and nodes[i + 1] not in pts:
                        pts.append(nodes[i + 1])
                        changed = True
            else:
                s = float(spacings.get(a, 1.0))
                t = float(tol.get(a, s / 4))
                if v <= lo + t:
                    nv = _clamp(a, lo - s)
                    if nv < lo - 1e-9:
                        pts.append(nv)
                        changed = True
                elif v >= hi - t:
                    nv = _clamp(a, hi + s)
                    if nv > hi + 1e-9:
                        pts.append(nv)
                        changed = True
                else:
                    left = max(p for p in pts if p <= v)
                    right = min(p for p in pts if p >= v)
                    if (right - left) > 2 * t and all(abs(v - p) > 1e-9 for p in pts):
                        pts.append(v)
                        changed = True
            pts.sort()

        if not changed:
            break

    return best, fitted, grid


def fit_all_lines(grid, lines, ref_grid=None, ref='LTE', integrate=True, limit=1.0, crop_kwargs=None):
    """Fit all lines in a grid and return the guesses."""
    crop_kwargs = crop_kwargs or {}
    ref_grid = ref_grid or grid
    sign = 1 if ref == 'LTE' else -1
    guesses = {}
    refs = {}
    for l0 in lines:
        line = grid.spectra[0].get_line(l0, verbose=False)
        kw = dict(crop_kwargs, xmin=line.lb, xmax=line.lr)
        refs[l0] = ref_grid.crop(da=0.0, LTE=(ref == 'LTE'), **kw)
        pfit, _ = grid.fit(refs[l0], fit=['da'], LTE=(ref != 'LTE'), integrate=integrate, **kw)
        guesses[l0] = np.clip(sign * pfit['da'], -limit, limit)
    return guesses, refs


def get_NLTE_correction(top, star, atom_file, atmos_file, lines, absmet_files=None, settings=None,
                        tol=0.01, limit=1.0, ref='LTE', integrate=True, verbose=False,
                        max_threads=None, launch_kwargs=None, crop_kwargs=None):
    """
    Compute NLTE corrections for a single star.

    Returns
    -------
    guesses : dict  {line_center: correction}
    runs    : dict  {da: Run}
    """
    launch_kwargs = launch_kwargs or {}
    crop_kwargs = crop_kwargs or {}

    if max_threads is None:
        max_threads = os.cpu_count() or 1

    atom_file = Path(atom_file)
    if not atom_file.exists():
        raise FileNotFoundError(f"Atom file {atom_file} not found!")

    sign = 1 if ref == 'LTE' else -1
    max_iter = 10

    teff, logg, feh = star[:3]
    vmic = star[3] if len(star) == 4 else 0.0

    name = atom_file.stem if atom_file.stem != 'atom' else atom_file.suffix[1:]
    name += f'_t{teff:04.0f}g{logg*10:02.0f}m{feh*10:+03.0f}'
    if vmic > 0:
        name += f'v{vmic:.0f}'

    if absmet_files is not None and feh in absmet_files:
        absmet_file = absmet_files[feh]
    else:
        print(f'No absmet file for [Fe/H]={feh}!')
        absmet_file = ''

    nml = {
        'timer_params': {'sec_per_report': 1e8, 'verbose': 2},
        'atmos_params': {'atmos_file': atmos_file, 'feh': feh, 'vmic': vmic},
        'atom_params': {
            'atom_file': atom_file,
            'relative_abundance': True,
            'use_atom_abnd': False,
        },
        'rt_iteration_params': {'maxiter': 99},
        'rt_output_params': {'maxiter': 3},
        'nlte_params': {'convlim': 0.01, 'convmax': 0.001},
        'opac_table_params': {'absmet_file': absmet_file},
    }
    if settings is not None:
        nml = update_dict(nml, settings)

    def build_grid():
        if any(runs[d] is None for d in runs):
            raise RuntimeError("Some runs failed, cannot build grid.")
        das = sorted(runs.keys())
        return tm.grid([runs[d] for d in das], da=das), das

    def find_interval(da_val, knots):
        idx = np.searchsorted(knots, da_val, side='right') - 1
        return max(0, min(idx, len(knots) - 2))

    # Step 1: Build initial coarse grid — all 3 runs in parallel
    da_grid = sorted([-limit, 0.0, limit])
    runs = _run_batch(_make_da_jobs(da_grid, nml, name, ref, top), launch_kwargs, max_threads)
    grid, _ = build_grid()

    refs = {}
    for l0 in lines:
        line = grid.spectra[0].get_line(l0, verbose=verbose)
        kw = dict(crop_kwargs, xmin=line.lb, xmax=line.lr)
        refs[l0] = grid.crop(da=0.0, LTE=(ref == 'LTE'), **kw)

    guesses, _ = fit_all_lines(grid, lines, ref=ref, crop_kwargs=crop_kwargs, integrate=integrate, limit=limit)

    for iteration in range(max_iter):
        prev_guesses = dict(guesses)
        knots = np.array(sorted(runs.keys()))

        intervals = defaultdict(list)
        for l0, g in guesses.items():
            ii = find_interval(g, knots)
            intervals[(ii, knots[ii], knots[ii + 1])].append((l0, g))

        new_points = []
        for (ii, lo, hi), line_guesses in intervals.items():
            if hi - lo < tol:
                continue
            guess_values = [g for _, g in line_guesses]
            split_at = np.median(guess_values)

            # Snap to nearest multiple of tol
            split_at = round(split_at / tol) * tol

            if (split_at - lo) < tol or (hi - split_at) < tol:
                split_at = (lo + hi) / 2.0
            if (split_at - lo) < tol or (hi - split_at) < tol:
                continue

            # Avoid exact duplicates of existing grid points
            if split_at not in runs:
                new_points.append(split_at)

        if not new_points:
            break

        # Run all new points in parallel
        new_runs = _run_batch(_make_da_jobs(new_points, nml, name, ref, top), launch_kwargs, max_threads)
        runs.update(new_runs)

        grid, _ = build_grid()
        guesses, _ = fit_all_lines(grid, lines, ref=ref, crop_kwargs=crop_kwargs, integrate=integrate, limit=limit)

        max_shift = max(abs(guesses[l0] - prev_guesses[l0]) for l0 in lines)

        if verbose:
            print(f"Iter {iteration}: added {len(new_points)} points, "
                  f"grid size={len(runs)}, max shift={max_shift:.2e}")

        if max_shift < tol:
            break

    return guesses, runs


def get_NLTE_corrections(top, stars, atom_file, atmos_files, lines, **kwargs):
    """ Compute NLTE corrections for multiple stars and return a DataFrame.

    Returns
    -------
    df : DataFrame with columns ['Teff', 'logg', '[Fe/H]', 'vmic', line1, line2, ...]
    runs : dict {star: {da: Run}} with the runs for each star.
    """
    atom_file = Path(atom_file)
    if not atom_file.exists():
        raise FileNotFoundError(f"Atom file {atom_file} not found!")

    df = pd.DataFrame()
    runs = {}

    for i, star in enumerate(stars):
        df.loc[i, ['Teff', 'logg', '[Fe/H]', 'vmic'][:len(star)]] = star

        if star not in atmos_files:
            print(f'No atmosphere file for star {star}!')
            continue
        else:
            atmos_file = atmos_files[star]

        corr, runs[star] = get_NLTE_correction(top, star, atom_file, atmos_file, lines, **kwargs)

        for l0 in lines:
            df.loc[i, l0] = corr[l0]

        if i == 0:
            print(df.to_string(float_format="{:.2f}".format))
        else:
            print('\n'.join(df.loc[i:i+1].to_string(float_format="{:.2f}".format).split('\n')[1:]))

    return df, runs
