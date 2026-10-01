"""Headless synthesis workflows shared by Python callers and front ends.

Wires the app to ``tumult.compute_spectrum``: given a target (Teff, logg, FeH, vmic) it
finds the enclosing model atmospheres in a configured grid, synthesises them in parallel
(``build_stellar_grid`` → shared dispatch service, disk-cached) and interpolates a
continuum-normalised ``Spectrum`` at the target — ready to overlay on the canvas.

The split of inputs (see the plan):
- **preferences** (``store``): the reusable namelist setup — grid folder/pattern/format,
  atom/abund/absmet/linelist files, solver/output/runtime knobs;
- **run settings** (right-sidebar panel): only Teff, logg, FeH, vmic, the abundances
  override string, the relative-abundances flag, and LTE/NLTE;
- every entry point takes a :class:`~tumult.synthesis.SynthRequest` and the run values an
  operation varies around it (:func:`tumult.synthesis.operations.run_values`): the request
  carries its wavelength windows (:func:`resolve_request`, :func:`execute_request`), while
  the batch functions (``compute*``) take one ``xmin/xmax`` pair or a list of windows.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import json
import os
import re
import shutil
from pathlib import Path

from .api import (compute_spectrum, compute_spectra,
                  compute_window_batch as _compute_window_batch,
                  normalise_windows, write_window_file,
                  PRESETS_DIR, preset, PROVENANCE_NAME)
from ..abundance import (
    SYNTH_PAD_AA, abund_to_xfe, convert_abundances,
    read_solar_abundances, xfe_to_abund,
)
from .cancel import CancelToken, CancelledError
from .inputs import Atmosphere, AtmosphereGrid, Snapshot
from .request import ResolvedSynthRequest, SynthesisResult, SynthRequest
from .snapshot_grid import SNAPSHOT_FORMATS, discover_grid
from ..auxfuncs import update_dict

# namelist presets now live with the core (tumult/presets) so the library and this app
# share one location; see tumult.synthesis.PRESETS_DIR.
_PRESETS = PRESETS_DIR

__all__ = ["discover_atmos_grid", "compile_pattern", "recover_params",
           "build_settings", "compute", "compute_batch", "compute_window_batch",
           "resolve_request", "execute_request",
           "solve_request_departures",
           "compute_reference",
           "read_reference_linelist",
           "abund_elements", "abund_table", "DEFAULT_PATTERN", "STAGGER_PATTERN",
           "FORMAT_PATTERNS", "default_pattern",
           "CancelToken", "CancelledError"]

# preset layered on top of base.nml for the 1D LTE reference: LTE (maxiter=0), atom-only
# (empty line list), relative abundances — see tumult/presets/corr_ref.nml.
_CORR_REF_PRESET = preset("corr_ref.nml")

# (key, UI label, input step, default) — the four interpolated stellar axes,
# shared by the Synthesize run form (synthpanel) and the Fit tab (fittab).
STELLAR_PARAMS = [("teff", "Teff [K]", 1, 5777), ("logg", "log g", 0.01, 4.4),
                  ("feh", "[Fe/H]", 0.01, 0.0), ("vmic", "Microturbulence ξ [km/s]", 0.1, 1.0)]

# the MARCS naming: p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_....mod
DEFAULT_PATTERN = "{Teff:d}_g{logg:+.1f}_m{:3}_t{vmic:02d}_st_z{FeH:+.2f}"
# the post-processed Stagger grid: t50g45m00/atm3d.t50g45m0005_00008_80x80x240 is
# Teff=5000, log g=4.5, [Fe/H]=0.0, model version 05, snapshot 8; the solar models write
# Teff in full (t5777g44m00), hence the first alternative. Two-digit names truncate the
# quarter-step log g values (g47 reads as 4.7). See the grid's INSTRUCTIONS.md.
STAGGER_PATTERN = ("atm3d.t{Teff:04d}g{logg*10:02d}{FeH*10:+03d}{}_{snap:05d}"
                   "|atm3d.t{Teff/100:02d}g{logg*10:02d}{FeH*10:+03d}{}_{snap:05d}")
# the M3DIS server layout: PLATO_t45.48g39.64m0.000_v1.0/m3dis_319.atmos is Teff=4548,
# log g=3.964, [Fe/H]=0.0, snapshot 319 (the reader finds the mesh beside it)
M3DIS_PATTERN = "_t{Teff/100:.2f}g{logg*10:.2f}m{FeH:.3f}_{}/m3dis_{snap:d}.atmos"

# the pattern a registered grid starts from, by atmos_format (case-insensitive)
FORMAT_PATTERNS = {"marcs": DEFAULT_PATTERN, "stagger2": STAGGER_PATTERN,
                   "m3dis": M3DIS_PATTERN}


def default_pattern(atmos_format):
    """The stock path pattern for ``atmos_format``, or ``""`` when none is known."""
    return FORMAT_PATTERNS.get(str(atmos_format or "").lower(), "")

# the stellar fields a pattern may bind; Teff/logg/FeH are required, vmic optional
_PATTERN_FIELDS = ("Teff", "logg", "FeH", "vmic", "snap")
_FIELD_NAMES = {name.lower(): name for name in _PATTERN_FIELDS}   # names are case-insensitive
_REQUIRED_FIELDS = ("Teff", "logg", "FeH")
# a signed integer/float/exponent token (what one stellar field captures in a filename)
_NUM_RE = r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"

# fallback resolving power when none is configured — high enough for line synthesis
_DEFAULT_R = 300000.0

# memoised grid scans: (folder, pattern, snapshots) -> (tree_mtime, grid); the MARCS folder
# has ~15k files, so re-globbing+parsing it on every synthesis is needless work.
_GRID_MEMO = {}
# folder -> (its mtime, its immediate subfolders), for _tree_mtime
_SUBFOLDERS = {}


# one {...} token in a pattern: an optional field name, an optional linear op (* / + -)
# with a number, and an optional :spec (format specs are ignored for matching; see
# compile_pattern for the matching keywords).
_TOKEN_RE = re.compile(
    r"\{(?P<name>[A-Za-z_]\w*)?"
    r"(?:\s*(?P<op>[*/+-])\s*(?P<num>\d+(?:\.\d+)?))?"
    r"(?::(?P<spec>[^}]*))?\}"
)
# top-level alternatives: a | outside any {...} token
_ALTERNATIVE_RE = re.compile(r"\|(?![^{]*\})")


def _transform(op, num):
    """Map a pattern's *physical → token* op to the inverse ``(scale, offset)`` so that
    ``physical = token * scale + offset`` recovers the physical value from a matched token.

    The pattern formats the physical value, so e.g. ``{Teff/100}`` writes ``token = Teff/100``
    ⇒ ``Teff = token * 100``. A bare field is the identity.
    """
    if op is None:
        return (1.0, 0.0)
    k = float(num)
    if op == "/":
        return (k, 0.0)
    if op == "*":
        return (1.0 / k, 0.0)
    if op == "+":
        return (1.0, -k)
    if op == "-":
        return (1.0, k)
    return (1.0, 0.0)


_FORMAT_RE = re.compile(
    r"^(?P<sign>[+ -])?(?P<zero>0)?(?P<width>\d+)?(?:\.(?P<precision>\d+))?(?P<type>[dfeEgG])?$")
_SIGN_RE = "[+\\-pPmM]"   # a sign as written (+/-) or as a name letter (p/m)


def _parse_spec(spec):
    """A token's ``:spec``, a Python format item, as its parsed parts (or ``None``)."""
    spec = (spec or "").strip()
    if not spec:
        return None
    fmt = _FORMAT_RE.match(spec)
    if fmt is None:
        raise ValueError(f"Unrecognised format item {spec!r} in pattern token")
    return fmt.groupdict()


class _Field:
    """How one bound token turns its matched text into a physical value.

    The format item is read as the one that *wrote* the token: a ``+`` flag makes the sign
    mandatory, a width fixes the digit count (less the sign when mandatory, less the
    decimals for ``f``), a precision fixes the decimals, ``d`` and ``.0f`` forbid them.
    Without a type the token is any signed number. ``p``/``m`` are accepted as ``+``/``-``.
    """

    def __init__(self, op, num, fmt):
        self.scale, self.offset = _transform(op, num)
        fmt = fmt or {}
        self.sign_required = fmt.get("sign") == "+"
        self.digits = None
        self.decimals = None
        kind = fmt.get("type")
        width = int(fmt["width"]) if fmt.get("width") else None
        if kind == "d":
            self.decimals = 0
        elif kind == "f":
            self.decimals = int(fmt["precision"]) if fmt.get("precision") is not None else None
        if width and kind in ("d", "f"):
            digits = width - (1 if self.sign_required else 0)
            if self.decimals:
                digits -= self.decimals + 1
            self.digits = digits if digits > 0 else None

    def regex(self, name):
        sign = _SIGN_RE + ("" if self.sign_required else "?")
        digits = rf"\d{{{self.digits}}}" if self.digits else r"\d+"
        if self.decimals is None:
            number = digits + r"(?:\.\d+)?(?:[eE][+-]?\d+)?"
        elif self.decimals == 0:
            number = digits
        else:
            number = digits + rf"\.\d{{{self.decimals}}}"
        return f"(?P<{name}>{sign}{number})"

    def value(self, raw):
        sign = 1.0
        if raw and raw[0] in "+-pPmM":
            sign = -1.0 if raw[0] in "-mM" else 1.0
            raw = raw[1:]
        # rounded so that inverting a scaled token (t65.26 -> 6526) lands on the value
        # the name was written from, which keys the node
        return round(sign * (float(raw) * self.scale + self.offset), 9) + 0.0


class _Pattern:
    """A compiled filename pattern: one regex per alternative, each with a named group per
    stellar field, plus the per-field rule recovering physical values from matched tokens."""

    def __init__(self, alternatives, fields):
        self.alternatives = alternatives  # [(regex, {field: _Field})]
        self.fields = fields              # ordered tuple of bound field names

    @property
    def has_vmic(self):
        return "vmic" in self.fields

    @property
    def has_snap(self):
        return "snap" in self.fields

    def search(self, name):
        """Recover ``{field: physical_value}`` from ``name``, or ``None`` if it doesn't match
        or a required field is missing/unparsable."""
        for regex, rules in self.alternatives:
            m = regex.search(name)
            if not m:
                continue
            out = {}
            for field, rule in rules.items():
                raw = m.group(field)
                if raw is None:
                    continue
                try:
                    out[field] = rule.value(raw)
                except ValueError:
                    return None
            if any(f not in out for f in _REQUIRED_FIELDS):
                return None
            return out
        return None


# compiled-pattern memo (patterns are tiny; this just avoids recompiling per file scan)
_PATTERN_MEMO = {}


def _compile_alternative(pattern):
    parts = []
    rules = {}
    pos = 0
    for m in _TOKEN_RE.finditer(pattern):
        parts.append(re.escape(pattern[pos:m.start()]))
        pos = m.end()
        name, op, num, spec = m.group("name", "op", "num", "spec")
        spec = _parse_spec(spec)
        name = _FIELD_NAMES.get((name or "").lower(), name)
        if name in _PATTERN_FIELDS:
            if name in rules:
                raise ValueError(f"Pattern binds {name!r} more than once: {pattern!r}")
            rules[name] = _Field(op, num, spec)
            parts.append(rules[name].regex(name))
        else:
            # unnamed / unknown field => skip; honour a fixed width from the spec
            width = (spec or {}).get("width")
            parts.append(f".{{{int(width)}}}" if width else ".*?")
    parts.append(re.escape(pattern[pos:]))
    missing = [f for f in _REQUIRED_FIELDS if f not in rules]
    if missing:
        raise ValueError(
            f"Pattern {pattern!r} must bind {', '.join(_REQUIRED_FIELDS)} "
            f"(missing: {', '.join(missing)}).")
    return re.compile("".join(parts)), rules


def compile_pattern(pattern):
    """Compile a filename ``pattern`` into a :class:`_Pattern`.

    A pattern reads like the Python format that wrote the name. Tokens are
    ``{Name[op NUM][:spec]}`` where ``Name`` is one of ``Teff``/``logg``/``FeH``/``vmic``/
    ``snap`` (case-insensitive; Teff/logg/FeH required, vmic and the snapshot number
    optional), ``op`` is one of ``* / + -`` applying a linear transform, and ``:spec`` is a
    Python format item. The format item constrains the match (see :class:`_Field`):
    ``{FeH*10:+03d}`` is a mandatory sign and two digits, ``{logg:+.1f}`` a sign and one
    decimal, ``{Teff:d}`` any integer. Signs may
    be written ``+``/``-`` or ``p``/``m``. An unnamed token (``{}`` or ``{:3}``) is a skip
    wildcard — a width consumes exactly that many characters, else lazily. Literal text
    is matched verbatim, anywhere in the file's path relative to the grid folder (head and
    tail are implicitly wild). ``|`` separates alternatives tried in order; every
    alternative must bind the same fields.

    So ``t{Teff/100:d}`` recovers Teff=5000 from ``t50``, and ``g{logg*10:02d}{FeH*10:+03d}``
    recovers (4.5, -3.0) from ``g45m30``. See :data:`STAGGER_PATTERN` for a complete example.
    """
    cached = _PATTERN_MEMO.get(pattern)
    if cached is not None:
        return cached
    alternatives = [_compile_alternative(text)
                    for text in _ALTERNATIVE_RE.split(pattern) if text.strip()]
    if not alternatives:
        raise ValueError("Empty pattern.")
    fields = tuple(alternatives[0][1])
    for _regex, rules in alternatives[1:]:
        if set(rules) != set(fields):
            raise ValueError(f"Pattern alternatives bind different fields: {pattern!r}")
    compiled = _Pattern(alternatives, fields)
    _PATTERN_MEMO[pattern] = compiled
    return compiled


def recover_params(pattern, name):
    """Recover ``(teff, logg, feh, vmic)`` from a single filename ``name`` using ``pattern``.

    Returns a 4-tuple (``vmic`` is ``None`` when the pattern has no vmic field), or ``None``
    if the name doesn't match. Used to pre-fill the run panel when a file is picked directly.
    """
    d = compile_pattern(pattern).search(name)
    if d is None:
        return None
    return (d["Teff"], d["logg"], d["FeH"], d.get("vmic"))


def _tree_mtime(folder):
    """The newest change of the folder or one of its immediate subfolders, which is what a
    recursive scan of a directory-per-node grid depends on.

    Adding or removing a subfolder changes the folder's own mtime, so the subfolder list
    is only listed again then; a flat grid of thousands of files costs a single stat.
    """
    own = folder.stat().st_mtime
    hit = _SUBFOLDERS.get(str(folder))
    if hit is None or hit[0] != own:
        # os.scandir reports each entry's type from the listing itself; Path.is_dir
        # would stat every file
        with os.scandir(folder) as entries:
            subfolders = [Path(entry.path) for entry in entries if entry.is_dir()]
        _SUBFOLDERS[str(folder)] = hit = (own, subfolders)
    return max([own, *(sub.stat().st_mtime for sub in hit[1])])


def _select_snapshots(items, policy, has_snap):
    """The files representing one node of a snapshot grid.

    ``items`` are ``(parent, snap, path)`` for every matching file of the node. When
    several subfolders map to the same node (an updated model beside its predecessor),
    the lexically last folder wins, as the last matching file wins for a 1D grid. Within
    it, ``policy`` (see :func:`snapshot_grid.parse_snapshot_selection`) picks by the
    bound snapshot number, or by name order when the pattern binds none; a count takes
    the leading entries of :func:`snapshot_grid.halton_order` over that sequence.
    """
    from .snapshot_grid import parse_snapshot_selection, select_count
    kind, arg = parse_snapshot_selection(policy)
    parent = max(item[0] for item in items)
    files = sorted(((snap, path) for folder, snap, path in items if folder == parent),
                   key=lambda item: (item[0] if item[0] is not None else 0, item[1].name))
    if kind == "all":
        chosen = files
    elif kind == "count":
        chosen = select_count(files, arg)
    else:
        if not has_snap:
            raise ValueError("Selecting snapshots by number needs a pattern that binds "
                             "{snap}.")
        chosen = [item for item in files if int(item[0]) in arg]
        missing = arg - {int(item[0]) for item in chosen}
        if missing:
            raise ValueError(f"{parent} has no snapshots {sorted(missing)}; available: "
                             f"{[int(item[0]) for item in files]}")
    return tuple(path for _snap, path in chosen)


def discover_atmos_grid(folder, pattern=DEFAULT_PATTERN, snapshots=None):
    """Scan ``folder`` for atmosphere files whose paths match ``pattern``.

    Files are listed recursively and matched on their path relative to ``folder``, so a
    flat MARCS folder matches on file names while a directory-per-node grid can bind its
    parameters from the folder name and the snapshot number from the file name.

    Returns ``{(teff, logg, feh, vmic): Path}`` when the pattern binds ``vmic``, else
    ``{(teff, logg, feh): Path}``; the last matching file (in sorted path order) wins a
    key. With ``snapshots`` (a selection policy, see :mod:`.snapshot_grid`) every key maps
    instead to the tuple of snapshot files representing that node. Cached per
    ``(folder, pattern, snapshots)`` and re-scanned when the folder tree's mtime changes.
    """
    folder = Path(folder).expanduser()
    if not folder.is_dir():
        raise FileNotFoundError(f"Atmosphere grid folder not found: {folder}")

    key = (str(folder), pattern, snapshots)
    mtime = _tree_mtime(folder)
    hit = _GRID_MEMO.get(key)
    if hit and hit[0] == mtime:
        return hit[1]

    p = compile_pattern(pattern)
    matches = {}
    for f in sorted(folder.rglob("*")):
        if not f.is_file():
            continue
        d = p.search(f.relative_to(folder).as_posix())
        if d is None:
            continue
        key2 = ((d["Teff"], d["logg"], d["FeH"], d["vmic"]) if p.has_vmic
                else (d["Teff"], d["logg"], d["FeH"]))
        matches.setdefault(key2, []).append((f.parent, d.get("snap"), f))
    if not matches:
        raise ValueError(f"No atmosphere files in {folder} matched pattern {pattern!r}")
    if snapshots is None:
        grid = {k: items[-1][2] for k, items in matches.items()}
    else:
        grid = {k: _select_snapshots(items, snapshots, p.has_snap)
                for k, items in matches.items()}
    _GRID_MEMO[key] = (mtime, grid)
    return grid


def _linelist_settings(path):
    """A ``linelist_params`` dict for a single shared file-or-folder input: the explicit
    files either way (:func:`tumult.linelist.linelist_files`), so dispatch never lists
    the folder itself."""
    if not path:
        return {}
    from ..linelist import linelist_files
    return {"linelist_file": [str(f) for f in linelist_files(path)]}


def build_settings(prefs, run_vals, atmos_format=None):
    """Assemble the nested nml ``settings`` dict (keyed by namelist section).

    ``prefs`` is the persisted preferences dict; ``run_vals`` holds the per-run values
    from the run panel (``teff/logg/feh/vmic/abundances/relative_abundances/mode`` and, for
    NLTE, ``atoms`` — a list of model-atom file paths — and optionally ``line_grid``, a
    generated bound-bound wavelength grid). ``atmos_format`` is the selected
    grid's format. The atmosphere file, FeH and vmic are filled per node by
    ``build_stellar_grid``.

    A 3D snapshot format (:data:`~tumult.synthesis.snapshot_grid.SNAPSHOT_FORMATS`)
    samples the background and the LTE line opacity from tables in both passes
    (``sample_back_opac`` and ``sample_linelist`` in ``rt_iteration_params`` and
    ``rt_output_params``). The flags go together: the Fortran allocates and builds the
    line table only in a pass that also samples the background, and a pass with
    ``sample_linelist`` alone reads no line opacity at all (``rt_patch_mod``). Departure
    species are masked out of the table and computed in detail, so NLTE lines are
    unaffected; a pass without a line list gets an empty table, which is how the departure
    phase of a snapshot runs. The caller merges ``config.expert`` on top of this dict, so an
    explicit value of either flag wins. The Fortran refuses an abundance ladder (a
    ``set_file`` with more than the unperturbed sets) under ``sample_linelist``, so a 3D
    ladder runs one process per rung (:func:`tumult.synthesis.api._foldable_ladder`).
    """
    nlte = str(run_vals.get("mode", "LTE")).upper() == "NLTE"
    # ``only_atom``: synthesise from the model atom's own transitions with NO external line
    # list (used by the abundance-correction tool so an LTE reference and an NLTE spectrum
    # differ only by LTE-vs-NLTE). ``atom_lines`` keeps the line list and adds the atom's
    # transitions to it (the Fortran drops the atom's species from the list). Either emits
    # the atom block below even in LTE mode.
    only_atom = bool(run_vals.get("only_atom", False))
    atom_lines = only_atom or bool(run_vals.get("atom_lines", False))

    atmos = {}
    fmt = atmos_format or prefs.get("atmos_format")
    if fmt:
        atmos["atmos_format"] = fmt
    if prefs.get("nx"): atmos["nx"] = int(prefs["nx"])
    if prefs.get("ny"): atmos["ny"] = int(prefs["ny"])
    if prefs.get("nz"): atmos["nz"] = int(prefs["nz"])

    composition = {}
    if prefs.get("abund_file"):
        composition["abund_file"] = prefs["abund_file"]
    if run_vals.get("abundances"):
        composition["abundances"] = run_vals["abundances"]
        composition["relative_abundances"] = bool(run_vals.get("relative_abundances", False))

    # atom-only runs carry no external line list (the corr_ref.nml preset also blanks it);
    # everything else uses the configured synthesis line list plus the element switches.
    linelist = {"linelist_file": [""], "linelist_folder": ""} if only_atom else _linelist_settings(prefs.get("linelist"))
    # element line-selection switches (Fortran /linelist_params/): a space-separated list of
    # element symbols. Used by the abundance pipeline's blend metric — synthesise a window
    # with the target element dropped (``exclude_elements``) or in isolation (``only_elements``).
    for key in ("exclude_elements", "only_elements"):
        val = run_vals.get(key)
        if val:
            linelist[key] = val if isinstance(val, str) else " ".join(val)

    # timer_params and the LTE/NLTE m3d_params/atom_params constants live in the
    # presets (tumult/presets/{base,nlte}.nml), applied as launch defaults.
    settings = {
        "linelist_params": linelist,
    }
    if str(fmt or "").lower() in SNAPSHOT_FORMATS:
        sampled = {"sample_back_opac": True, "sample_linelist": True}
        settings["rt_iteration_params"] = dict(sampled)
        settings["rt_output_params"] = dict(sampled)
    # A curve-of-growth rung only repeats one atmosphere at another abundance; the
    # transition dump is the same diagnostic every time, so a caller synthesising a
    # ladder asks for it once. Absent means the Fortran default (write it).
    if run_vals.get("save_linelist") is False:
        settings.setdefault("rt_output_params", {})["save_linelist"] = False
    if atmos:
        settings["atmos_params"] = atmos
    if composition:
        settings["composition_params"] = composition

    if nlte or atom_lines:
        atoms = run_vals.get("atoms") or []
        if isinstance(atoms, str):
            atoms = [atoms]
        # a one-element list still serialises to a Fortran array (see NameList); the
        # Fortran atom_params namelist accepts an array of atom_file entries.
        atom = {"atom_file": [str(a) for a in atoms]}
        # a generated bound-bound grid (tumult.line_quadrature) replaces the atoms' own
        # quadratures; the atom copy that goes with it expresses its windows in QNORM_KMS
        if run_vals.get("line_grid"):
            from ..line_quadrature import QNORM_KMS
            atom["line_grid_file"] = str(run_vals["line_grid"])
            atom["qnorm"] = float(QNORM_KMS)
        if prefs.get("convlim"):
            atom["convlim"] = float(prefs["convlim"])
        if prefs.get("convmax"):
            atom["convmax"] = float(prefs["convmax"])
        settings["atom_params"] = atom

    return settings


def _absmet_files(prefs, atmos_files):
    """Map every grid [Fe/H] to the configured absmet file (None if unset).

    A single ``absmet`` file applies to all metallicities; a folder is parsed per-FeH.
    Returning a per-FeH dict is required because ``build_stellar_grid`` overrides
    ``composition_params.absmet_file`` per node (and would otherwise blank it).
    """
    absmet = prefs.get("absmet")
    if not absmet:
        return None
    p = Path(absmet).expanduser()
    fehs = {round(k[2], 4) for k in atmos_files}
    if p.is_dir():
        found = {}
        for k in atmos_files:
            feh = k[2]
            cand = sorted(p.glob(f"*{feh:+.2f}*"))
            if cand:
                found[feh] = str(cand[0])
        return found or None
    return {k[2]: str(p) for k in atmos_files} if fehs else None


def _prune_spec(prefs, prune):
    """The per-node line-list pruning descriptor for the API, or None for no pruning.

    Splits the decision in two: the skirt is a preference (``linelist_skirt_aa``, 0 to switch
    pruning off), because how far outside a window a line may still matter is a property of
    the physics and not of any one caller; the composition the selection is made at is the
    caller's (see :func:`tumult.synthesis.api._node_linelists`). A caller passing
    ``prune=None`` gets no pruning whatever the preference says.
    """
    if prune is None:
        return None
    skirt = float(prefs.get("linelist_skirt_aa") or 0.0)
    if not skirt > 0:
        return None
    return {"skirt": skirt, "abundances": str(prune.get("abundances") or ""),
            "relative_abundances": bool(prune.get("relative_abundances", True))}


def _check_run_values(values):
    """Reject the run values no execution path supports, naming the cause."""
    if values.get("only_atom") and not values.get("atoms"):
        raise ValueError("Atom-only synthesis needs model atoms.")
    if str(values.get("mode", "LTE")).upper() == "NLTE" and not values.get("atoms"):
        raise ValueError("NLTE synthesis needs model atoms.")


def _stellar(request):
    """The request's star as the ``(teff, logg, feh, vmic)`` node key; a snapshot without
    microturbulence carries 0."""
    star = request.star
    return (star.teff, star.logg, star.feh, star.vmic if star.vmic is not None else 0.0)


def _prepare(request, run_values_list, extra_presets=None, *, set_spectrum=True,
             window_batch=False):
    """Resolve common inputs once; each backend retains its own window assembly."""
    from .operations import run_values
    config, atmosphere = request.config, request.atmosphere
    prefs, root = config.to_prefs(), config.root
    star = _stellar(request)
    if isinstance(atmosphere, AtmosphereGrid):
        grid = atmosphere.record(config)
        atmos_files = discover_grid(grid, root)
        mesh = _grid_mesh(grid)
    else:
        atmos_files = {star[:3]: config.resolve(atmosphere.path)}
        mesh = {}
    values_list = [run_values(request, **dict(extra)) for extra in run_values_list]
    for values in values_list:
        _check_run_values(values)
    modes = {str(values.get("mode", "LTE")).upper() for values in values_list}
    if len(modes) != 1:
        raise ValueError("A synthesis batch must share one LTE/NLTE mode.")
    nlte = modes.pop() == "NLTE"
    R = config.resolving_power
    settings_list = []
    for values in values_list:
        settings = build_settings(prefs, values, atmos_format=atmosphere.format)
        settings = update_dict(settings, {"atmos_params": mesh})
        settings = update_dict(settings, deepcopy(config.expert))
        if set_spectrum:
            sampling = settings.setdefault("spectrum_params", {})
            if window_batch or not any(sampling.get(k) for k in ("daa", "nnu", "lam_file")):
                sampling["R_sampling"] = R
        settings_list.append(settings)
    defaults = [str(_PRESETS / "base.nml")]
    if nlte:
        defaults.append(str(_PRESETS / "nlte.nml"))
    defaults.extend(str(path) for path in (extra_presets or []))
    launch_kwargs = {"overwrite": False, "verbose": 0, "default": defaults,
                     "cwd": str(root)}
    # This value only bounds this batch's own outstanding jobs (below, and in
    # _run_batch's per-batch allocation); it must not mutate the shared execution
    # service's global budget, which is resolved once by the daemon at startup
    # (plan/multi-processing.md, "CPU budget").
    max_threads = config.threads or (os.cpu_count() or 1)
    common = dict(top=str(config.cache), star=star,
                  atmos_files=atmos_files, absmet_files=_absmet_files(prefs, atmos_files),
                  method=config.interpolation,
                  max_threads=max_threads, launch_kwargs=launch_kwargs,
                  vmic_step=config.vmic_step, nlte=nlte)
    return common, settings_list, R


def compute(request, run_values=None, xmin=None, xmax=None, *, extra_presets=None,
            set_spectrum=True, cancel_token=None, return_runs=False):
    """Compute the normalised spectrum of ``request`` with ``run_values`` over ``[xmin, xmax]``.

    An atmosphere grid is interpolated at the request's star; an individual atmosphere is
    synthesised directly (a one-node grid keyed at the star, weight 1.0). The request's own
    windows are not read.

    ``extra_presets`` (optional) is a list of namelist preset paths appended after the
    base/NLTE presets, so later ones win — used to layer ``corr_ref.nml`` for the 1D LTE
    reference synthesis.

    ``set_spectrum`` (default True) injects a sampling trigger into ``spectrum_params`` (the
    configured R) so the Fortran emits a spectrum. The 1D LTE reference passes it False so the
    sampling is left entirely to the ``corr_ref.nml`` preset.

    Blocking — call from a worker thread. ``build_stellar_grid`` fans synthesis of the
    enclosing atmospheres out across processes internally. ``cancel_token`` (optional) lets
    a caller hard-kill the dispatch runs; a cancel raises ``CancelledError``.
    """
    common, settings, _ = _prepare(request, [run_values or {}], extra_presets,
                                    set_spectrum=set_spectrum)
    return compute_spectrum(**common, settings=settings[0], xmin=xmin, xmax=xmax,
                            cancel_token=cancel_token, return_runs=return_runs)


def compute_batch(request, run_values_list, xmin=None, xmax=None, *, extra_presets=None,
                  set_spectrum=True, cancel_token=None, fold=False, store=None,
                  abundances=None):
    """Batched :func:`compute`: synthesise several run-value variants of ``request`` (e.g. a
    curve-of-growth's abundance ladder), fanning *every* dispatch run through the shared
    execution service via :func:`tumult.synthesis.compute_spectra`. Returns a list of
    Spectrum aligned with ``run_values_list`` (empty list for empty input).

    Each variant is assembled through the same :func:`build_settings` + preset + window path as
    :func:`compute`, so a variant's runs are byte-identical to the ``compute`` call it replaces
    — the batch only changes scheduling (all points run concurrently instead of one at a time),
    not results. All variants must share the synthesis mode (LTE/NLTE) and preset/line-list
    config; only per-point knobs (typically ``abundances``) should differ. Every variant
    keeps its own run and departure solve unless an LTE caller opts into ``fold=True`` (see
    :func:`tumult.synthesis.compute_spectra`). ``store`` and ``abundances`` keep the runs in
    memory for a workflow (see :class:`tumult.synthesis.api.NodeStore`).
    """
    if not run_values_list:
        return []
    common, settings, _ = _prepare(request, run_values_list, extra_presets,
                                    set_spectrum=set_spectrum)
    return compute_spectra(**common, settings_list=settings, xmin=xmin, xmax=xmax,
                           cancel_token=cancel_token, fold=fold, store=store,
                           abundances=abundances)


def node_selector(request, store):
    """``select(params) -> frozenset`` of the synthesis node keys a star on ``request``'s
    atmosphere grid interpolates from, selected as :func:`compute_window_batch` selects
    them and kept in ``store`` (a :class:`~tumult.synthesis.api.NodeStore`).

    The spatial selection does not depend on vmic and the vmic lattice not on the other
    parameters, so each is kept separately: stars that differ only in vmic share one
    spatial selection.
    """
    from .atmos_select import _vmic_levels
    config = request.config
    atmos_files = discover_grid(request.atmosphere.record(config), config.root)
    absmet_files = _absmet_files(config.to_prefs(), atmos_files)
    spatial = {}

    def select(params):
        position = tuple(float(params[k]) for k in ("teff", "logg", "feh"))
        vmic = float(params["vmic"])
        if position not in spatial:
            nodes, _ = store.select((*position, vmic), atmos_files,
                                    method=config.interpolation,
                                    vmic_step=config.vmic_step, absmet_files=absmet_files,
                                    vmic_points=2)
            # a snapshot node has no vmic lattice; its key stands as it is
            spatial[position] = (frozenset(n[0][0] for n in nodes if n[4] is not None),
                                 frozenset(n[0] for n in nodes if n[4] is None))
        lattice, fixed = spatial[position]
        levels = [v for v, _w in _vmic_levels(vmic, config.vmic_step)]
        return fixed | {(site, v) for site in lattice for v in levels}
    return select


def compute_window_batch(request, run_values_list, windows, *, extra_presets=None,
                         cancel_token=None, prune=None, fold=False, store=None,
                         abundances=None, probes=()):
    """Batched :func:`compute` over an abundance ladder *and* several windows at once.

    Same contract as :func:`compute_batch` for the ladder axis, with ``windows`` a
    sequence of ``(xmin, xmax)`` bounds instead of one window. The ladder folds onto the
    run's output patches; each variant is its own run unless ``fold=True`` puts an LTE
    ladder onto one run's set axis. Returns ``results[variant][window]``.

    Falls back to one synthesis per window when the variants cannot fold, so callers can
    use this wherever they would otherwise loop :func:`compute_batch` over windows.

    ``prune`` asks for the per-node line-list pruning: the ``{"abundances",
    "relative_abundances"}`` the selection is made at (see :func:`_prune_spec`).

    Unfolded variants are still scheduled together, and a run keyed on one abundance is
    reused by any later ladder that asks for the same lattice rung. ``store`` and
    ``abundances`` keep those runs in memory for a workflow (see
    :class:`tumult.synthesis.api.NodeStore`). ``probes`` are further stellar parameters
    (``Star`` or mappings) interpolated from the same runs; the result then has one entry
    per star, the request's first.
    """
    if not run_values_list or not windows:
        return []
    common, settings, R = _prepare(request, run_values_list, extra_presets,
                                    window_batch=True)
    return _compute_window_batch(**common, settings_list=settings, windows=windows,
                                 R=R, cancel_token=cancel_token,
                                 prune=_prune_spec(request.config.to_prefs(), prune), fold=fold,
                                 store=store, abundances=abundances,
                                 probes=[_stellar(request.update(star=probe))
                                         for probe in probes])


def _grid_mesh(grid):
    """The ``nx/ny/nz`` a registered snapshot grid asks its cubes to be sampled onto, when
    it names one; otherwise the preferences' global mesh stands."""
    resolution = (grid or {}).get("resolution")
    if not resolution:
        return {}
    values = [int(v) for v in resolution]
    if len(values) != 3 or any(v <= 0 for v in values):
        raise ValueError("A grid resolution needs three positive integers (nx, ny, nz).")
    return dict(zip(("nx", "ny", "nz"), values))


# --- resolved public request -----------------------------------------------------------
def _resolved_path(root, value):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def _preparation_grid(request):
    """The MARCS grid record HOP preparations draw their reference atmospheres from: the
    atmosphere's ``preparation`` grid, else a MARCS grid itself, else none."""
    atmosphere, config = request.atmosphere, request.config
    if atmosphere.preparation is not None:
        return atmosphere.preparation.record(config)
    if isinstance(atmosphere, AtmosphereGrid) and atmosphere.is_marcs:
        return atmosphere.record(config)
    return {}


def _nearest_grid_path(atmospheres, stellar, vmic=0.0):
    keys = sorted({tuple(float(v) for v in key[:3]) for key in atmospheres})
    spans = [max(values) - min(values) or 1.0 for values in zip(*keys)]
    target = tuple(float(v) for v in stellar[:3])
    nearest = min(keys, key=lambda key: (sum(((a - b) / span) ** 2
                                             for a, b, span in zip(key, target, spans)), key))
    candidates = [(key, path) for key, path in atmospheres.items()
                  if tuple(float(v) for v in key[:3]) == nearest]
    return min(candidates,
               key=lambda item: (abs((float(item[0][3]) if len(item[0]) > 3 else 0.0)
                                    - float(vmic)), str(item[1])))[1]


def nearest_marcs_atmosphere(grid, target, prefs):
    """The MARCS model of ``grid`` nearest to ``target`` (teff, logg, feh[, vmic]).

    Used by preparations that need one representative 1D model of a star, such as the
    bound-bound grid for its NLTE corrections; ``grid`` must be a MARCS-format grid record.
    """
    if str(grid.get('format', 'Marcs')).lower() != 'marcs':
        raise ValueError('A MARCS-format grid is needed to pick a reference model.')
    root = execution_root(prefs).resolve()
    atmospheres = discover_atmos_grid(_resolved_path(root, grid['folder']),
                                      grid.get('pattern') or DEFAULT_PATTERN)
    vmic = float(target[3]) if len(target) > 3 else 0.0
    return str(_resolved_path(root, _nearest_grid_path(atmospheres, target, vmic)))


def preparation_root(prefs):
    """Where HOP preparations (pruned line lists, bound-bound grids) are cached."""
    root = execution_root(prefs).resolve()
    return _resolved_path(root, prefs.get('synth_output_dir') or 'synth') / '_preparation'


def _reference_atmospheres(request, prefs, root, prefs_vmic=1.0):
    """The MARCS atmospheres HOP preparations use, and which synthesis atmosphere each serves.

    An individual MARCS atmosphere and the nodes of a MARCS grid are prepared on
    themselves. Other formats use the nearest node of their preparation grid.
    """
    atmosphere, star = request.atmosphere, request.star
    if isinstance(atmosphere, Atmosphere) and atmosphere.is_marcs:
        source = str(_resolved_path(root, atmosphere.path))
        return (source,), {source: source}
    grid = _preparation_grid(request)
    if not grid.get('folder'):
        return (), {}
    if str(grid.get('format', 'Marcs')).lower() != 'marcs':
        raise ValueError('The automatic HOP preparation reference grid must use MARCS format.')
    references = discover_atmos_grid(
        _resolved_path(root, grid['folder']), grid.get('pattern') or DEFAULT_PATTERN)
    if isinstance(atmosphere, AtmosphereGrid):
        from .atmos_select import select_nodes
        synthesis_grid = discover_grid(atmosphere.record(request.config), root)
        stellar = star.tuple() if star.vmic is not None else star.tuple()[:3]
        nodes, _ = select_nodes(stellar, synthesis_grid,
                                method=prefs.get('interp_method', 'NaturalNeighbor'),
                                vmic_step=float(prefs.get('vmic_step') or 0.5))
        aliases = {}
        for node in nodes:
            source = str(_resolved_path(root, node[6]))
            # a snapshot node has no vmic; its MARCS reference takes the request's
            vmic = node[4] if node[4] is not None else float(prefs_vmic)
            reference = str(_resolved_path(
                root, _nearest_grid_path(references, node[1:4], vmic)))
            aliases[source] = reference
        return tuple(dict.fromkeys(aliases.values())), aliases

    vmic = star.vmic if star.vmic is not None else 0.0
    reference = str(_resolved_path(root, _nearest_grid_path(references, star.tuple()[:3],
                                                            vmic)))
    source = str(_resolved_path(root, atmosphere.path))
    return (reference,), {source: reference}


def _single_atmosphere_settings(request, prefs):
    """What a grid node's runs are given per node (:func:`tumult.synthesis.api._node_atmos`,
    the absmet file of its [Fe/H]) for an individual atmosphere or snapshot, whose runs
    bypass the node selection when departures are solved or read back directly."""
    atmosphere, star = request.atmosphere, request.star
    if isinstance(atmosphere, AtmosphereGrid):
        return {}
    feh = float(star.feh)
    settings = {}
    absmet = _absmet_files(prefs, {(star.teff, star.logg, feh): atmosphere.path}) or {}
    if absmet.get(feh):
        settings['opac_table_params'] = {'absmet_file': absmet[feh]}
    if isinstance(atmosphere, Atmosphere):
        settings['atmos_params'] = {'feh': feh, 'vmic': float(star.vmic)}
    return settings


def resolve_request(request, run_values=None) -> ResolvedSynthRequest:
    """Resolve ``request`` with ``run_values`` (see
    :func:`tumult.synthesis.operations.run_values`): paths, windows, namelist settings and
    preparation references. Raises ``ValueError`` for run values no execution path
    supports."""
    if not isinstance(request, SynthRequest):
        raise TypeError('resolve_request expects a SynthRequest')
    from .operations import run_values as complete_values
    # a JSON round-trip owns the inputs: nothing stays shared with the caller
    request = SynthRequest.from_dict(json.loads(request.to_json()))
    extra = json.loads(json.dumps(dict(run_values or {})))
    values = complete_values(request, **extra)
    _check_run_values(values)
    config, atmosphere = request.config, request.atmosphere
    root, prefs = config.root, config.to_prefs()
    settings = build_settings(prefs, values, atmos_format=atmosphere.format)
    if isinstance(atmosphere, AtmosphereGrid):
        settings = update_dict(settings,
                               {'atmos_params': _grid_mesh(atmosphere.record(config))})
    settings = update_dict(settings, _single_atmosphere_settings(request, prefs))
    settings = update_dict(settings, deepcopy(config.expert))
    # A correction ladder without windows uses the model atom's frequency grid,
    # including a prepared line_grid_file. Injecting R alone requests an explicit
    # spectrum grid without supplying its bounds, which DISPATCH cannot construct.
    atom_grid = not request.windows and (values.get('atom_lines') or values.get('only_atom'))
    if not atom_grid:
        settings.setdefault('spectrum_params', {}).setdefault('R_sampling', config.resolving_power)
    if settings.get('spectrum_params', {}).get('R_sampling', _DEFAULT_R) <= 0:
        raise ValueError('Resolving power must be positive.')
    windows = request.windows
    padded = tuple((window[0] - SYNTH_PAD_AA, window[1] + SYNTH_PAD_AA)
                   for window in windows)
    references, aliases = _reference_atmospheres(
        request, prefs, root, prefs_vmic=float(values.get('vmic') or 1.0))
    if (not values.get('only_atom') and config.linelist_skirt > 0 and config.linelist
            and not references):
        raise ValueError('Automatic line-list preparation has no MARCS reference grid. '
                         'Configure a reference grid or set linelist_skirt=0 to use the full list.')
    return ResolvedSynthRequest(request, extra, settings, windows, padded, references,
                                aliases, str(root), str(output_root(prefs)))


def _prepare_artifact(resolved, cancel_token=None):
    values = resolved.values
    if values.get("only_atom"):
        return None
    config, windows = resolved.request.config, resolved.padded_windows
    prefs = config.to_prefs()
    skirt = config.linelist_skirt
    if not (skirt > 0 and windows and resolved.preparation_references and prefs.get('linelist')):
        return None
    from ..linestrength import prepare_linelist
    settings = resolved.settings.get('composition_params', {}) or {}
    return prepare_linelist(preparation_root(prefs),
                            resolved.preparation_references, prefs['linelist'], windows,
                            R=float(resolved.settings['spectrum_params']['R_sampling']), dlam=skirt,
                            vmic=float(values.get('vmic', 1.0)),
                            feh=float(values.get('feh', 0.0)),
                            abundances=settings.get('abundances', ''),
                            relative_abundances=bool(settings.get('relative_abundances', True)),
                            aliases=resolved.preparation_aliases,
                            cwd=resolved.execution_root, cancel_token=cancel_token)


def execute_request(request, run_values=None, *, cancel_token=None):
    """Synthesise ``request`` with ``run_values``, or a :class:`ResolvedSynthRequest`, for
    every atmosphere kind through one result contract."""
    from ..departures import DepartureGrid, Departures
    resolved = (request if isinstance(request, ResolvedSynthRequest)
                else resolve_request(request, run_values))
    request, atmosphere = resolved.request, resolved.request.atmosphere
    config, values = request.config, resolved.values
    settings = deepcopy(resolved.settings)
    artifact = _prepare_artifact(resolved, cancel_token)
    grid = isinstance(atmosphere, AtmosphereGrid)
    if artifact is not None and not grid:
        settings.setdefault('linelist_params', {}).pop('linelist_folder', None)
        settings['linelist_params']['linelist_file'] = [artifact.files[0]]
    launch_kwargs = {'overwrite': False, 'verbose': 0,
                     'default': [str(_PRESETS / 'base.nml')], 'cwd': resolved.execution_root}
    nlte = str(values.get('mode', 'LTE')).upper() == 'NLTE'
    if nlte:
        launch_kwargs['default'].append(str(_PRESETS / 'nlte.nml'))
    max_threads = config.threads or (os.cpu_count() or 1)
    top = str(config.cache)
    path = None if grid else _resolved_path(Path(resolved.execution_root), atmosphere.path)

    spectra, runs = [], []
    # requested departures are solved inside the runs themselves; only solved ones are read
    departures = (request.departures
                  if isinstance(request.departures, (Departures, DepartureGrid)) else None)
    dep_dirs = departures.folders() if isinstance(departures, DepartureGrid) else None
    if isinstance(atmosphere, Snapshot) or (departures is not None and dep_dirs is None):
        from .snapshot import solve_departures, synthesise, synthesise_snapshot
        resolution = atmosphere.resolution if isinstance(atmosphere, Snapshot) else None
        if nlte and departures is None:
            departures = solve_departures(top, path, settings, resolution=resolution,
                                          atmos_format=atmosphere.format,
                                          launch_kwargs=launch_kwargs, max_threads=max_threads,
                                          cancel_token=cancel_token)
        # one spectrum run covers every window through a window file, as on the 1D path;
        # the departure solve above never sees it, since its grid comes from the model atom
        item = deepcopy(settings)
        if resolved.windows:
            window_file = write_window_file(
                Path(resolved.execution_root) / top / '_windows',
                normalise_windows(resolved.windows, SYNTH_PAD_AA))
            item.setdefault('spectrum_params', {})['window_file'] = str(window_file)
        if nlte:
            run = synthesise(departures, settings=item, top=top, launch_kwargs=launch_kwargs,
                             max_threads=max_threads, cancel_token=cancel_token)
        else:
            run = synthesise_snapshot(top, path, settings=item, resolution=resolution,
                                      atmos_format=atmosphere.format,
                                      launch_kwargs=launch_kwargs, max_threads=max_threads,
                                      cancel_token=cancel_token)
        runs.append(str(Path(run.folder).resolve()))
        spectra.extend(run.crop(window=window) if window is not None else run.crop()
                       for window in (resolved.windows or (None,)))
    else:
        extra = [resolved.run_values]
        if resolved.windows:
            common, _settings, R = _prepare(request, extra, window_batch=True)
            spectra_by_run, run_paths = _compute_window_batch(
                **common, settings_list=[settings], windows=resolved.windows, R=R,
                pad=SYNTH_PAD_AA, cancel_token=cancel_token,
                preparation=artifact if grid else None, return_runs=True, dep_dirs=dep_dirs)
            spectra.extend(spectra_by_run[0])
            runs.extend(run_paths)
        else:
            common, _settings, _ = _prepare(request, extra)
            value, run_paths = compute_spectrum(**common, settings=settings,
                                                cancel_token=cancel_token, return_runs=True,
                                                dep_dirs=dep_dirs)
            spectra.append(value)
            runs.extend(str(path) for path in run_paths)
    windows = resolved.windows or (None,)
    return SynthesisResult(tuple(spectra), windows, request, resolved.run_values, settings,
                           artifact, departures, tuple(runs))


def solve_request_departures(request, run_values, *, cancel_token=None):
    """Solve only the departures of ``request`` with NLTE ``run_values`` on its individual
    atmosphere or snapshot; the request's windows and line-list pruning do not enter."""
    config = replace(request.config, linelist_skirt=0.0)
    resolved = resolve_request(request.update(config=config, windows=()), run_values)
    request, atmosphere = resolved.request, resolved.request.atmosphere
    if isinstance(atmosphere, AtmosphereGrid):
        raise ValueError('Choose an individual 1D atmosphere or 3D snapshot for saved departures.')
    if str(resolved.values.get('mode', 'LTE')).upper() != 'NLTE':
        raise ValueError('Departure solving requires NLTE mode.')
    defaults = [str(_PRESETS / 'base.nml'), str(_PRESETS / 'nlte.nml')]
    launch_kwargs = {'overwrite': False, 'verbose': 0, 'default': defaults,
                     'cwd': resolved.execution_root}
    from .snapshot import solve_departures
    resolution = atmosphere.resolution if isinstance(atmosphere, Snapshot) else None
    return solve_departures(str(request.config.cache),
                            _resolved_path(Path(resolved.execution_root), atmosphere.path),
                            resolved.settings, resolution=resolution,
                            atmos_format=atmosphere.format, launch_kwargs=launch_kwargs,
                            max_threads=request.config.threads or (os.cpu_count() or 1),
                            cancel_token=cancel_token)


# --- fitting ---------------------------------------------------------------------------
def fitted_abundance_terms(best, elements, *, base="", relative, feh, solar):
    """The abundance override string describing a fit result, in the caller's convention.

    ``best`` is the fit's parameter dict, whose abundance axes are always fitted as [X/Fe];
    ``relative`` says which convention the caller states abundances in, so an absolute
    caller gets A(X) back rather than the [X/Fe] the fit worked in. ``base`` is any
    fixed override string the fit held constant, kept ahead of the fitted terms.
    """
    terms = [str(base).strip()] if str(base).strip() else []
    for el in elements:
        if el in best:
            xfe = float(best[el])
            val = xfe if relative else xfe_to_abund(xfe, el, feh, solar)
            terms.append(f"{el}={val:+.2f}")
    return " ".join(terms)


def abund_table(prefs):
    """``{element: A(X)}`` from the configured abundance file (solar reference).

    Thin prefs adapter: resolves the (possibly ``dispatch_cwd``-relative) path from the
    app preferences and delegates the parsing to
    :func:`tumult.abundance.read_solar_abundances`. Returns ``{}`` if the file is unset
    or unreadable.
    """
    path = prefs.get("abund_file")
    if not path:
        return {}
    p = Path(path).expanduser()
    if not p.is_absolute() and prefs.get("dispatch_cwd"):
        p = Path(str(prefs["dispatch_cwd"])).expanduser() / p
    return read_solar_abundances(p)


def abund_elements(prefs):
    """Element symbols available for the "+" dropdown (excludes H and Fe, the ``[Fe/H]`` axis)."""
    return [el for el in abund_table(prefs) if el not in ("H", "Fe")]


# --- 1D LTE reference (abundance corrections) ------------------------------------------

def compute_reference(request, atom, xmin=None, xmax=None, *, cancel_token=None,
                      return_runs=False, line_grid=None, with_linelist=False):
    """Synthesise the 1D LTE reference spectrum of a model ``atom`` for ``request``.

    The request's atmosphere (a grid interpolated at its star, or one atmosphere) and
    composition are synthesised in LTE and atom-only via the ``corr_ref.nml`` preset
    (``maxiter=0``, empty line list); the request's departures do not enter. The atom's own
    bound-bound transitions are the lines, and the run also writes ``linelist_out.txt``,
    read back by :func:`read_reference_linelist`. ``line_grid`` (optional) is a generated
    bound-bound grid the run samples instead of the atom's own quadrature. ``with_linelist``
    adds the configured synthesis line list (blends from the other species) to the atom's
    transitions. Returns the normalised reference :class:`Spectrum`. Blocking; call from a
    worker thread.
    """
    run_values = {"only_atom": not with_linelist, "atom_lines": True, "atoms": [str(atom)]}
    if line_grid:
        run_values["line_grid"] = str(line_grid)
    return compute(request.update(departures=None), run_values, xmin, xmax,
                   extra_presets=[_CORR_REF_PRESET], set_spectrum=False,
                   cancel_token=cancel_token, return_runs=return_runs)


def read_reference_linelist(run_dirs):
    """Read the atom transitions from the runs that produced a reference spectrum.

    Each node synthesises the same atom transitions; the first producing node is
    representative. Missing or malformed output is an error, never an unrelated cache hit.
    """
    from ..linelist import LineList

    if not run_dirs:
        raise ValueError("The reference spectrum has no producing runs.")
    dump = Path(run_dirs[0]) / "linelist_out.txt"
    frame = LineList.from_tumult_format(str(dump)).reset_index()
    return frame[[c for c in ("ion", "lam0", "evi", "loggf", "strength", "i", "j")
                  if c in frame.columns]]


# --- run cache management --------------------------------------------------------------
# dispatch writes node runs under `data/<synth_output_dir>/...` relative to the app cwd;
# each completed run is a folder containing params.nml (see compute_spectrum for the path
# scheme: <output_dir>/<window>_<settings-hash>/t..g..m..[v..]).
_RUN_RE = re.compile(r"^t(\d+)g(\d+)m([+-]\d+)(?:v(\d+))?$")


def execution_root(prefs):
    """Directory against which dispatch resolves its relative inputs and outputs."""
    return Path(prefs.get("dispatch_cwd") or ".").expanduser().resolve()


def output_root(prefs):
    """Directory the synthesis outputs and the ``_workflows/`` records live under."""
    return _resolved_path(execution_root(prefs), prefs.get("synth_output_dir") or "synth")


def cache_root(prefs):
    """Directory holding cached MULTI3D node runs."""
    return execution_root(prefs) / "data" / (prefs.get("synth_output_dir") or "synth")


def input_root(prefs):
    """Directory holding the input namelists launch() writes (parallel to the outputs).

    launch() writes each job's params namelist at ``<synth_output_dir>/<bucket>/<run>.nml``
    relative to the app cwd, while the run outputs land under ``data/<synth_output_dir>/``
    (see cache_root). Clearing the cache must remove both trees.
    """
    return execution_root(prefs) / (prefs.get("synth_output_dir") or "synth")


def _parse_run(name):
    """Approximate stellar params from a run folder name (lossy: names are rounded).

    Only for a run with no ``provenance.json`` -- one launched before that file existed.
    """
    m = _RUN_RE.match(name)
    if not m:
        return None
    teff, glog, mh, vmic = m.groups()
    # vmic is encoded in tenths of km/s (e.g. 'v15' -> 1.5 km/s)
    return {"teff": int(teff), "logg": int(glog) / 10,
            "feh": int(mh) / 10, "vmic": int(vmic) / 10 if vmic is not None else 0}


def _provenance(run_dir):
    """``(star, abund)`` a run recorded for itself, or ``(None, None)``.

    :func:`tumult.synthesis.api._save_provenance` writes this beside every run it launches,
    holding the values that built the namelist. A run from before it existed has none, and
    then the folder name is all there is -- which is rounded, and carries no abundances.
    """
    try:
        record = json.loads((Path(run_dir) / PROVENANCE_NAME).read_text())
    except (OSError, ValueError):
        return None, None
    abund = record.get("abund") or {}
    return record.get("star"), {"items": [(el, float(v)) for el, v in abund.get("items", [])],
                                "relative": bool(abund.get("relative", False))}


def _window(bucket):
    """(lo, hi) Å from a 'w<lo>_<hi>[...]' bucket/run name, else None."""
    m = re.match(r"^w(\d+)_(\d+)", bucket or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


def _dir_size(path):
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def list_cache(prefs):
    """All cached runs, most recent first. Each: path, run, star, abund, kind, window, size, mtime.

    ``kind`` is ``"departure"`` for NLTE departure-coefficient runs (the ``dep_<hash>`` bucket,
    no output spectrum) and ``"spectrum"`` for the synthesised spectra.
    """
    root = cache_root(prefs)
    if not root.is_dir():
        return []
    entries = []
    for nml in root.rglob("params.nml"):
        run_dir = nml.parent
        rel = run_dir.relative_to(root).parts
        bucket = rel[0] if len(rel) > 1 else ""
        run = rel[-1]
        star, abund = _provenance(run_dir)
        entries.append({
            "path": str(run_dir.resolve()),
            "run": run,
            "star": star or _parse_run(run),   # the run's own record, else its rounded name
            "abund": abund or {"items": [], "relative": False},
            "kind": "departure" if bucket.startswith("dep") else "spectrum",
            "window": _window(bucket) or _window(run),
            "size": _dir_size(run_dir),
            "mtime": run_dir.stat().st_mtime,
        })
    entries.sort(key=lambda e: e["mtime"], reverse=True)
    return entries


def delete_cache_entry(prefs, path):
    """Delete one cached run folder (guarded to stay inside the cache root)."""
    root = cache_root(prefs).resolve()
    p = Path(path).resolve()
    if root not in p.parents or not (p / "params.nml").is_file():
        return False
    shutil.rmtree(p)
    bucket = p.parent  # drop the now-empty window/hash bucket too
    if bucket != root and bucket.is_dir() and not any(bucket.iterdir()):
        bucket.rmdir()

    # remove the matching input namelist (mirrors the run path under input_root)
    in_root = input_root(prefs).resolve()
    rel = p.relative_to(root)
    in_nml = (in_root / rel).with_suffix(".nml")
    in_nml.unlink(missing_ok=True)
    in_bucket = in_nml.parent  # drop the now-empty bucket too
    if in_bucket != in_root and in_bucket.is_dir() and not any(in_bucket.iterdir()):
        in_bucket.rmdir()
    return True


def clear_cache(prefs):
    """Delete every cached run under the cache root. Returns the number removed."""
    root = cache_root(prefs).resolve()
    if not root.is_dir():
        return 0
    runs = list_cache(prefs)
    for child in list(root.iterdir()):
        shutil.rmtree(child) if child.is_dir() else child.unlink()

    # also clear the parallel input-namelist tree (see input_root)
    in_root = input_root(prefs)
    if in_root.resolve() != root and in_root.is_dir():
        for child in list(in_root.iterdir()):
            shutil.rmtree(child) if child.is_dir() else child.unlink()
    return len(runs)
