"""LineMasks — the shared fit/synth-window model for TumultStudio2.

A :class:`LineMask` is one wavelength window used to fit synthesis to a target. It carries a
*fit window* ``[lmin, lmax]`` (the χ²/EW region — for imported masks, exactly the file's
base/top columns) and a wider *synth/visual window* ``[smin, smax]`` (the RT synthesis span
so wings compute, the visual context, and the region a local continuum is fit over). One mask
may be the home of several bound-bound transitions, so ``lam0`` is a single representative
centre (``NaN`` when the file gives none).

A :class:`LineMaskSet` is an ordered list of masks plus the ``pad`` used to derive each mask's
synth window from its fit window. It is the single currency both the abundance pipeline and the
Synthesize/Fit tab consume; it is produced either by auto-generation (master line list + target)
or by importing a linemask file.

Linemask files (see ``input/linemasks/``) come in several flavours — this module reads them all:

* header rows (``wave_peak  wave_base  wave_top  note`` or ``lmin lmax [lam0]``) or none,
* ``;`` / ``#`` comment lines, and inline trailing ``; comment`` after data,
* 2 columns (``lmin lmax`` segments), 3 (``lam0 lmin lmax``) or 4 (extra species code / note),
* tab / comma / whitespace separators, and nm-or-Å wavelengths (auto-converted to Å).

The canonical *save* format is ``lam0  lmin  lmax  [note]`` (peak-first).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, fields
from pathlib import Path

import astropy.units as u
import numpy as np


# Å added each side of the fit window to form the synth/visual window on import or edit;
# a LineMaskSet may override it. Distinct from the abundance pipeline's SYNTH_PAD_AA, which
# pads *this* window again before synthesis: a line is computed over
# ``fit_window ± (DEFAULT_PAD + SYNTH_PAD_AA)``. The two currently share a value but answer
# different questions — how much context a window shows, and how much extra the wings need —
# so they are calibrated apart.
DEFAULT_PAD = 0.5

# below this wavelength (in the file's own units) an optical linemask must be in nm, not Å
_NM_THRESHOLD = 2000.0

_NAN = float("nan")

# header-name synonyms → canonical role
_LAM0_KEYS = {"lam0", "wave_peak", "peak", "wave_center", "center", "centre", "wave0"}
_LMIN_KEYS = {"lmin", "wave_base", "base", "lo", "low", "wmin", "left", "start"}
_LMAX_KEYS = {"lmax", "wave_top", "top", "hi", "high", "wmax", "right", "end"}
_NOTE_KEYS = {"note", "species", "element", "ion", "label", "comment"}

# a token counts as a wavelength/number only if it parses as a float
_FLOAT_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def _is_float(tok: str) -> bool:
    return bool(_FLOAT_RE.match(tok))


@dataclass
class LineMask:
    lam0: float = _NAN          # representative centre [Å] (NaN if none / multi-transition)
    lmin: float = _NAN          # fit-window lower bound [Å]
    lmax: float = _NAN          # fit-window upper bound [Å]
    smin: float = _NAN          # synth/visual-window lower bound [Å]
    smax: float = _NAN          # synth/visual-window upper bound [Å]
    note: str = ""              # species / label; may host several transitions
    include: bool = True        # count this mask when fitting / aggregating
    renorm: bool = False        # retained in old mask records; measurements do not apply it
    cont_pts: tuple = ()        # retained raw-flux markers from old mask records

    _TUPLE_FIELDS = ("cont_pts",)

    # --- window views (so callers can keep using fit_window / window) ------------------
    @property
    def fit_window(self) -> tuple:
        return (self.lmin, self.lmax)

    @property
    def synth_window(self) -> tuple:
        return (self.smin, self.smax)

    # alias: the abundance pipeline historically called the wide window ``window``
    @property
    def window(self) -> tuple:
        return (self.smin, self.smax)

    def pad_synth(self, pad: float) -> None:
        """Derive the synth/visual window as the fit window widened by ``pad`` Å each side."""
        self.smin = self.lmin - pad
        self.smax = self.lmax + pad

    def set_fit(self, lo: float, hi: float, pad: float) -> None:
        """Set ``[lmin, lmax]`` (ascending) and re-derive the synth window from ``pad``."""
        self.lmin, self.lmax = (float(lo), float(hi)) if lo <= hi else (float(hi), float(lo))
        self.pad_synth(pad)

    def renorm_knots(self) -> list | None:
        """The 4 colinear continuum knots ``[(x,y), …]`` this mask injects into the object's
        spline continuum when :attr:`renorm` is on, or ``None`` when it has no markers.

        The two inner markers (``cont_pts``) define a straight line; the two outer knots sit at
        ``smin`` / ``smax`` on that same line (linear extrapolation), so the whole ``[smin, smax]``
        continuum segment is perfectly linear while the padding beyond it blends into the rest of
        the spline. See :func:`extrapolate_edge_knots`.
        """
        if not (self.renorm and len(self.cont_pts) == 2):
            return None
        (x1, y1), (x2, y2) = self.cont_pts
        return extrapolate_edge_knots((x1, y1), (x2, y2), self.smin, self.smax)

    # --- serialisation -----------------------------------------------------------------
    def to_dict(self) -> dict:
        d = {f.name: getattr(self, f.name) for f in fields(self)}
        d["cont_pts"] = [list(p) for p in self.cont_pts]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "LineMask":
        names = {f.name for f in fields(cls)}
        d = {k: v for k, v in (d or {}).items() if k in names}
        if d.get("cont_pts"):
            d["cont_pts"] = tuple(tuple(float(c) for c in p) for p in d["cont_pts"])
        return cls(**d)

    def to_row(self) -> str:
        """One canonical ``lam0  lmin  lmax  [note]`` line (peak-first)."""
        lam0 = self.lam0 if np.isfinite(self.lam0) else 0.5 * (self.lmin + self.lmax)
        cells = [f"{lam0:10.4f}", f"{self.lmin:10.4f}", f"{self.lmax:10.4f}"]
        if self.note:
            cells.append(self.note)
        return "  ".join(cells)


@dataclass
class LineMaskSet:
    masks: list = field(default_factory=list)   # list[LineMask]
    pad: float = DEFAULT_PAD
    source_path: str = ""
    label: str = ""
    uid: int = 0                # identity in a multi-set MaskStore (0 = unassigned)
    # the element these windows diagnose ("Fe", "Eu", …), set by hand — guessed from the
    # file name on import. It is deliberately the *only* species information a set carries:
    # a window generally spans several transitions, so per-window ion / E_low / log gf are
    # not well defined here and are recovered downstream from a synthesis run's
    # ``linelist_out.txt`` (which weights each transition by its actual contribution).
    element: str = ""
    # pristine copy of the imported (or last-saved) state, so the GUI can tell an edited set
    # from the file it came from and offer a revert. Empty for a set drawn from scratch.
    baseline: dict = field(default_factory=dict, repr=False)

    @property
    def name(self) -> str:
        """Display name for the sidebar: the file-derived label, else a placeholder."""
        return self.label or "Untitled masks"

    # --- edited-since-import tracking --------------------------------------------------
    def _state(self) -> dict:
        """The comparable content of the set (its masks and padding); excludes identity
        and the baseline itself, so it can be diffed against :attr:`baseline`."""
        return {"masks": [m.to_dict() for m in self.masks],
                "pad": self.pad, "element": self.element}

    def mark_clean(self) -> None:
        """Adopt the current content as the pristine state (on import and after a save)."""
        self.baseline = self._state()

    @property
    def dirty(self) -> bool:
        """Has the set been edited away from the file it was imported from / saved to?"""
        return bool(self.baseline) and self._state() != self.baseline

    def revert(self) -> bool:
        """Restore the pristine state; False when there is nothing to revert to."""
        if not self.baseline:
            return False
        self.masks[:] = [LineMask.from_dict(d) for d in self.baseline["masks"]]
        self.pad = float(self.baseline["pad"])
        self.element = str(self.baseline.get("element") or "")
        return True

    def __iter__(self):
        return iter(self.masks)

    def __len__(self):
        return len(self.masks)

    def to_dict(self) -> dict:
        return {"masks": [m.to_dict() for m in self.masks], "pad": self.pad,
                "source_path": self.source_path, "label": self.label, "uid": self.uid,
                "element": self.element, "baseline": self.baseline}

    @classmethod
    def from_dict(cls, d: dict) -> "LineMaskSet":
        d = dict(d or {})
        return cls(masks=[LineMask.from_dict(m) for m in d.get("masks", [])],
                   pad=float(d.get("pad", DEFAULT_PAD)),
                   source_path=d.get("source_path", ""), label=d.get("label", ""),
                   uid=int(d.get("uid", 0)), element=str(d.get("element") or ""),
                   baseline=dict(d.get("baseline") or {}))

    def repad(self, pad: float) -> None:
        """Re-derive every mask's synth window from a new ``pad`` (Å)."""
        self.pad = float(pad)
        for m in self.masks:
            m.pad_synth(self.pad)


# --- element identification ----------------------------------------------------------------

# ``; element: Fe`` in a saved file's header; the only metadata we round-trip through the
# text format (the parser skips ``;`` lines, so older files simply have none)
_ELEMENT_RE = re.compile(r"^\s*[;#]\s*element\s*[:=]\s*([A-Za-z]{1,3})\s*$", re.I)

# names carry the element and little else: "fe-lmask.txt" -> Fe, "Fe/mask_test.txt" -> Fe.
# These are stripped before the leading token is tested against the periodic table.
_NAME_NOISE = re.compile(r"(?:^|[_\-.])(lmask|linemask|mask|masks|seg|segments|lines|line|"
                         r"test|edited|new|final|v\d+)(?=$|[_\-.])", re.I)


def is_element(symbol) -> bool:
    """Is ``symbol`` a chemical element symbol ("Fe" yes, "fe"/"Xx" no)?"""
    from tumult.isotopes import element_z
    try:
        element_z(str(symbol))
    except (KeyError, TypeError):
        return False
    return True


def guess_element(path) -> str:
    """Best guess at the element a line-mask file diagnoses, from its name then its folder.

    ``input/linemasks/Fe/fe-lmask.txt`` -> ``"Fe"``. Returns ``""`` when neither the file
    stem nor the parent directory names an element — the user then sets it by hand.
    """
    path = Path(path)
    for raw in (path.stem, path.parent.name):
        token = _NAME_NOISE.sub("", str(raw)).strip("_-. ")
        token = re.split(r"[_\-.\s]", token)[0] if token else ""
        if token:
            cand = token.capitalize()
            if is_element(cand):
                return cand
    return ""


def read_element_header(text: str) -> str:
    """The ``; element: X`` declaration from a saved file's comment header, or ``""``."""
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped[0] not in ";#":
            break                     # header is over — data has started
        m = _ELEMENT_RE.match(stripped)
        if m and is_element(m.group(1).capitalize()):
            return m.group(1).capitalize()
    return ""


# --- file parsing -------------------------------------------------------------------------

def _split_cells(line: str) -> list[str]:
    """Tokens of one data line: strip an inline ``;`` / ``#`` trailing comment, then split on
    commas, tabs or whitespace runs."""
    for c in (";", "#"):
        i = line.find(c)
        if i != -1:
            line = line[:i]
    line = line.strip()
    if not line:
        return []
    return re.split(r"[,\t]|\s+", line)


def _header_map(cells: list[str]) -> dict | None:
    """If ``cells`` is a header row (no leading float), map column index → role; else ``None``."""
    if not cells or _is_float(cells[0]):
        return None
    roles = {}
    for i, c in enumerate(cells):
        key = c.strip().lower()
        if key in _LAM0_KEYS:
            roles[i] = "lam0"
        elif key in _LMIN_KEYS:
            roles[i] = "lmin"
        elif key in _LMAX_KEYS:
            roles[i] = "lmax"
        elif key in _NOTE_KEYS:
            roles[i] = "note"
    return roles


def _parse_row(cells: list[str], roles: dict | None):
    """Parse one data row into ``(lam0, lmin, lmax, note)`` in the file's own units, or ``None``.

    With a ``roles`` header map, columns are read by name. Without one the leading float tokens
    are the wavelengths in the consistent peak-first convention of ``input/linemasks/`` (2 →
    ``lmin,lmax``; ≥3 → ``lam0, lmin, lmax`` positionally, so a centre offset from its window —
    e.g. an Hβ wing mask — is read correctly; any 4th float is a species code folded into the
    note) and the trailing non-float tokens are the note (e.g. ``Ni 1``).
    """
    if roles:
        num_idx = {i: r for i, r in roles.items() if r in ("lam0", "lmin", "lmax")}
        nums = {r: float(cells[i]) for i, r in num_idx.items()
                if i < len(cells) and _is_float(cells[i])}
        if "lmin" not in nums or "lmax" not in nums:
            return None
        # every column that is not a mapped wavelength is note text (keeps e.g. "Co 1" whole)
        note = " ".join(cells[i] for i in range(len(cells)) if i not in num_idx).strip()
        lam0 = nums.get("lam0", _NAN)
        lo, hi = nums["lmin"], nums["lmax"]
        return lam0, min(lo, hi), max(lo, hi), note

    # headerless: leading floats are wavelengths, the rest is the note
    floats, note_parts, in_note = [], [], False
    for tok in cells:
        if not in_note and _is_float(tok):
            floats.append(float(tok))
        else:
            in_note = True
            note_parts.append(tok)
    note = " ".join(note_parts).strip()
    if len(floats) >= 4:
        # lam0,lmin,lmax + a species code (e.g. 12.01) → keep the code in the note
        code = floats[3]
        note = (note + f" {code:g}").strip()
        floats = floats[:3]
    if len(floats) >= 3:
        lam0, b, c = floats[:3]                   # peak-first: col0 is the centre
        return lam0, min(b, c), max(b, c), note
    if len(floats) == 2:
        lo, hi = sorted(floats)
        return _NAN, lo, hi, note
    return None


def _detect_nm(values: list[float]) -> bool:
    """True when the wavelengths are plainly in nm (all well below the optical-Å floor)."""
    finite = [v for v in values if np.isfinite(v) and v > 0]
    return bool(finite) and float(np.median(finite)) < _NM_THRESHOLD


def parse_linemask_file(path, pad: float = DEFAULT_PAD, label: str = "") -> LineMaskSet:
    """Read a linemask file into a :class:`LineMaskSet` (windows in Å).

    Robust to the format variety in ``input/linemasks/``: comment lines, inline comments,
    optional headers, 2/3/4 columns, mixed separators, ``\\r`` line endings, and nm-vs-Å
    wavelengths (converted to Å when the file is plainly in nm). Each mask's synth/visual window
    is derived as the fit window ± ``pad`` Å.
    """
    path = Path(path)
    text = path.read_text()   # universal-newline translation folds \r / \r\n into \n
    roles, rows = None, []
    seen_data = False
    for raw_line in text.splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped[0] in ";#":
            continue                              # blank or full-line comment
        cells = _split_cells(raw_line)
        if not cells:
            continue
        if not seen_data:
            hm = _header_map(cells)
            if hm is not None:
                roles = hm
                seen_data = True                  # header consumed; only one allowed
                continue
        seen_data = True
        parsed = _parse_row(cells, roles)
        if parsed is not None:
            rows.append(parsed)

    # unit pass: convert nm → Å across the whole file when needed
    flat = [v for lam0, lo, hi, _ in rows for v in (lam0, lo, hi)]
    if _detect_nm(flat):
        rows = [(_conv(lam0), _conv(lo), _conv(hi), note) for lam0, lo, hi, note in rows]

    masks = []
    for lam0, lo, hi, note in rows:
        m = LineMask(lam0=lam0, lmin=lo, lmax=hi, note=note)
        m.pad_synth(pad)
        masks.append(m)
    # the element is declared in the header of files we wrote; otherwise guess from the name
    return LineMaskSet(masks=masks, pad=pad, source_path=str(path),
                       label=label or path.stem,
                       element=read_element_header(text) or guess_element(path))


def _conv(v: float) -> float:
    return float((v * u.nm).to_value(u.AA)) if np.isfinite(v) else v


def write_linemask_file(path, mset: LineMaskSet, included_only: bool = True) -> None:
    """Write ``mset`` in the canonical ``lam0  lmin  lmax  [note]`` format.

    Only the masks the user kept switched on are written (``included_only``): saving is how a
    generated candidate list becomes a clean, reusable mask file, so the per-mask toggle
    curates the output.
    """
    path = Path(path)
    masks = [m for m in mset.masks if m.include] if included_only else list(mset.masks)
    element = f"; element: {mset.element}\n" if mset.element else ""
    header = (f"; linemask — {mset.label}\n" + element
              + "; lam0        lmin        lmax      note\n")
    body = "\n".join(m.to_row() for m in masks)
    path.write_text(header + body + "\n")


# --- local-continuum helpers (Phase 6 renorm) ---------------------------------------------

def extrapolate_edge_knots(p1: tuple, p2: tuple, smin: float, smax: float) -> list:
    """The four colinear continuum knots for a mask: the two inner markers plus edge knots at
    ``smin`` / ``smax`` placed on the marker line (linear extrapolation). Returned x-sorted."""
    (x1, y1), (x2, y2) = p1, p2
    slope = (y2 - y1) / (x2 - x1) if x2 != x1 else 0.0
    y_lo = y1 + slope * (smin - x1)
    y_hi = y1 + slope * (smax - x1)
    knots = [(smin, y_lo), (x1, y1), (x2, y2), (smax, y_hi)]
    return sorted(knots, key=lambda p: p[0])


def seed_cont_markers(raw, mask: LineMask, pct: float = None) -> tuple:
    """Initial ``((x1,y1),(x2,y2))`` inner markers for a mask's local continuum.

    A straight line is fit to the *upper*-``pct``-percentile raw-flux points inside the visual
    window (the continuum, above the line core), and the two markers are that line evaluated at
    the fit-window edges — so the seeded continuum already follows a sloped continuum. Falls back
    to a flat percentile level when the window is too sparse.
    """
    from .continuum import SEED_PERCENTILE, flat_level, upper_envelope_fit

    pct = SEED_PERCENTILE if pct is None else pct
    xx, yy = raw.values
    xx = np.asarray(xx, float)
    yy = np.asarray(yy, float)
    inwin = (xx >= mask.smin) & (xx <= mask.smax)
    xi, yi = xx[inwin], yy[inwin]
    x1 = float(mask.lmin if np.isfinite(mask.lmin) else mask.smin)
    x2 = float(mask.lmax if np.isfinite(mask.lmax) else mask.smax)
    fit = upper_envelope_fit(xi, yi, pct)
    if fit is not None:
        slope, intercept = fit
        return ((x1, slope * x1 + intercept), (x2, slope * x2 + intercept))
    level = flat_level(yi, pct) if xi.size else 1.0
    return ((x1, level), (x2, level))
