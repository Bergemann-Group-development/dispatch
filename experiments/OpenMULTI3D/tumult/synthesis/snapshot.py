"""Solve departures on one 3D snapshot, then synthesise from them as often as wanted.

The stellar-grid path already works this way: :func:`tumult.synthesis.api.build_stellar_grid`
solves departures per atmosphere node with no line list and no formal solution, caches them
under a digest that deliberately excludes the window and the line list, and then runs the
spectrum step against the saved populations with ``maxiter=0``. Pointed at a single file it
will even accept a snapshot, since ``atmos_file`` and the destination ``nx``/``ny``/``nz``
both pass through it untouched. Three things make that the wrong shape for a cube:

* it crosses every atmosphere with a microturbulence lattice and writes ``vmic`` into each
  node's ``atmos_params``, which imposes a microturbulence on a snapshot that already
  carries a velocity field — and, off a lattice point, solves the same cube twice;
* it selects and interpolates between nodes bracketing a (teff, logg, feh) point, which one
  cube is not in;
* it keeps the solved departures to itself. ``dep_dirs`` is a local of one call, so a
  caller who wants a second window off the same solve has to reconstruct settings that
  digest identically and trust the cache, rather than holding the thing it solved.

This module is that same two-phase split with a snapshot for an input:

    field = solve_departures(top, snapshot, settings, resolution=(30, 30, 128))
    blue  = synthesise(field, window=(5160, 5190))
    red   = synthesise(field, window=(5520, 5540))       # no second NLTE solve

Both phases go through the same ``_run_batch`` the grid path uses, so they share the CPU
budget, the cancellation contract and the in-flight de-duplication with every other job in
the process. The on-disk format is unchanged: the departures phase writes the ``save``
directory (``atmos_mesh.txt``, ``pop_<atom>_meta.txt`` and its binaries) that
``read_external_mod``/``nlte_patch_mod::init_departures`` already read, so nothing in the
Fortran changes.

Coarsening the snapshot is the existing ``atmos_params`` ``nx``/``ny``/``nz`` point sampling
through ``read_atmos_mod::interpolate_atmos``; ``resolution`` is those three numbers and
nothing more. Users choose the sampling and assess its effect on the departures.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from tumult.auxfuncs import update_dict
from tumult.departures import Departures, atom_ids

from .api import _dep_bucket
from .executor import _run_batch

# What the departures phase must not carry: a departure solve has no output spectrum, so a
# line list or a wavelength window in its settings would both slow it down and split its
# cache by something it does not depend on. This is the whole point of the two-phase split.
SPECTRUM_ONLY = ('linelist_params', 'spectrum_params')


def _spectrum_name(nml):
    """Cache identity of a spectrum run: the whole namelist.

    Not ``_dep_bucket``'s digest, which deliberately ignores the line list and the window —
    the two things that distinguish one spectrum from another off the same departures.
    """
    digest = hashlib.sha1(
        json.dumps(nml, sort_keys=True, default=str).encode()).hexdigest()[:8]
    return f'spec_{digest}'


def departure_settings(settings, snapshot, resolution=None, atmos_format=None):
    """The namelist of the departures phase: the snapshot, no spectrum, saved populations."""
    nml = copy.deepcopy(settings or {})
    for section in SPECTRUM_ONLY:
        nml.pop(section, None)
    atmos = {'atmos_file': str(snapshot), 'save_atmos': True}
    if atmos_format:
        atmos['atmos_format'] = atmos_format
    if resolution is not None:
        atmos.update(zip(('nx', 'ny', 'nz'), (int(n) for n in resolution)))
    nml = update_dict(nml, {
        'atmos_params': atmos,
        'rt_output_params': {'angle_quad': 'none'},
    })
    # The populations are only worth saving from an NLTE solve; an expert's explicit
    # convlim and iteration limit are kept.
    nml.setdefault('nlte_params', {}).setdefault('convlim', 1e-3)
    return nml


def spectrum_settings(field, settings=None, window=None):
    """The namelist of the spectrum phase: the saved mesh and populations, no solve.

    It starts from the departure run's own settings rather than from nothing: frozen
    departures are departures *of* an atom, so the spectrum has to be told the same
    ``atom_params`` the solve used, and a caller who had to restate them could restate them
    wrongly. ``settings`` is merged on top — the line list, the sampling, an abundance
    override — and the phase's own three decisions are applied last, since they are what
    makes this the spectrum phase rather than another solve.
    """
    nml = copy.deepcopy(field.settings)
    # The destination mesh and the snapshot come back from the saved run's params.nml, so
    # restating them here would be a second and silently divergent statement of both.
    nml.pop('atmos_params', None)
    output = nml.get('rt_output_params')
    if output and str(output.get('angle_quad', '')).strip() == 'none':
        # the departures phase turned the formal solution off; this phase is the formal solution
        output.pop('angle_quad')
        if not output:
            nml.pop('rt_output_params')

    nml = update_dict(nml, copy.deepcopy(settings or {}))
    for axis in ('nx', 'ny', 'nz'):
        nml.get('atmos_params', {}).pop(axis, None)
    nml = update_dict(nml, {
        'atmos_params': {'atmos_format': 'save', 'atmos_file': str(field.folder),
                         'read_departures': True},
        'rt_iteration_params': {'maxiter': 0},
    })
    if window is not None:
        lo, hi = float(min(window)), float(max(window))
        nml = update_dict(nml, {'spectrum_params': {'aa_blue': lo, 'aa_red': hi}})
    return nml


def solve_departures(top, snapshot, settings=None, *, resolution=None, atmos_format=None,
                     launch_kwargs=None, max_threads=None, cancel_token=None):
    """Solve NLTE populations on ``snapshot`` and return the :class:`Departures`.

    The run lands in a ``dep_<digest>`` directory under ``top``, digested over the
    departure-relevant namelist sections only (:func:`~tumult.synthesis.api._dep_bucket`),
    so a second call with the same physics is served from the cache and a call that changes
    only the window or the line list never reaches this phase at all.
    """
    nml = departure_settings(settings, snapshot, resolution, atmos_format)
    name = _dep_bucket(top, nml)
    runs = _run_batch([('departures', name, nml)], launch_kwargs or {},
                      max_threads or 1, cancel_token=cancel_token)
    run = runs.get('departures')
    if run is None:
        raise RuntimeError(f"departure solve failed for {snapshot}; see {name}")
    folder = Path(run.folder)
    actual = tuple(int((nml.get('atmos_params') or {}).get(k))
                   for k in ('nx', 'ny', 'nz')
                   if (nml.get('atmos_params') or {}).get(k) is not None)
    field = Departures(folder=folder, name=Path(name), atmosphere=str(snapshot),
                           resolution=actual, species=atom_ids(folder), settings=nml,
                           composition=copy.deepcopy(nml.get('composition_params', {}) or {}))
    field.save()
    return field


def synthesise(field, *, settings=None, window=None, name=None, top=None, launch_kwargs=None,
               max_threads=None, cancel_token=None):
    """Synthesise one spectrum from an already-solved :class:`Departures`.

    Only the formal solution runs: the populations are read back and frozen. ``settings``
    carries whatever the spectrum needs and the departures did not — the line list, the
    sampling, an abundance override. An abundance override still scales the LTE populations
    the frozen coefficients multiply, so the line strength does follow it; whether the
    coefficients themselves may be left at another abundance's values is a scientific
    choice for the user to assess.

    The run goes under ``top``, the caller's cache, and beside the departures without one.
    """
    if not field.saved:
        raise ValueError(f"{field.folder} holds no saved departures to synthesise from")
    nml = spectrum_settings(field, settings, window)
    top = Path(top) if top is not None else field.name.parent
    # The spectrum's identity is its whole namelist, departure folder included, so a second
    # window or line list against the same field is a different run and a repeat is a hit.
    name = top / (name or _spectrum_name(nml))
    runs = _run_batch([('spectrum', name, nml)], launch_kwargs or {},
                      max_threads or 1, cancel_token=cancel_token)
    run = runs.get('spectrum')
    if run is None:
        raise RuntimeError(f"spectrum run failed for {field.folder}; see {name}")
    return run


def synthesise_snapshot(top, snapshot, *, settings=None, resolution=None,
                        atmos_format=None, name=None, launch_kwargs=None,
                        max_threads=None, cancel_token=None):
    """Run one snapshot directly, without stellar nodes or a vmic lattice."""
    nml = copy.deepcopy(settings or {})
    atmos = {'atmos_file': str(snapshot)}
    if atmos_format:
        atmos['atmos_format'] = atmos_format
    if resolution is not None:
        atmos.update(zip(('nx', 'ny', 'nz'), (int(v) for v in resolution)))
    nml = update_dict(nml, {'atmos_params': atmos})
    name = Path(top) / (name or _spectrum_name(nml))
    runs = _run_batch([('spectrum', name, nml)], launch_kwargs or {},
                      max_threads or 1, cancel_token=cancel_token)
    run = runs.get('spectrum')
    if run is None:
        raise RuntimeError(f'snapshot synthesis failed for {snapshot}; see {name}')
    return run
