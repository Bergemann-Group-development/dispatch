"""Parser + round-trip tests for the LineMask model (TumultStudio2/linemask.py).

Every file under ``input/linemasks/`` must parse into sane, ascending-bound masks with
wavelengths normalised to Å, and the canonical writer must round-trip them.
"""
import math
from pathlib import Path

import numpy as np
import pytest

from tumult import linemask
from tumult.linemask import LineMask, LineMaskSet, parse_linemask_file, write_linemask_file

ROOT = Path(__file__).resolve().parent.parent
LMASK_DIR = ROOT / "input" / "linemasks"
ALL_FILES = sorted(LMASK_DIR.rglob("*.txt"))


def test_have_linemask_files():
    assert ALL_FILES, f"no linemask files under {LMASK_DIR}"


@pytest.mark.parametrize("path", ALL_FILES, ids=lambda p: str(p.relative_to(LMASK_DIR)))
def test_parses_to_sane_masks(path):
    mset = parse_linemask_file(path)
    assert isinstance(mset, LineMaskSet)
    # a few files are entirely commented out (all data rows disabled with ';') — allow empty,
    # but every mask that IS produced must be well-formed and in Å.
    for m in mset.masks:
        assert np.isfinite(m.lmin) and np.isfinite(m.lmax)
        assert m.lmin < m.lmax
        assert m.smin < m.lmin and m.smax > m.lmax          # padded synth window encloses fit
        assert 1000.0 < m.lmin < 100000.0                   # optical/NIR Å, never left in nm
        # note: lam0 need NOT lie within [lmin,lmax] — e.g. Hβ wing masks sit off-centre.


def test_nm_file_converted_to_angstrom():
    # Ni values are ~481 nm in the file → ~4812 Å after conversion
    mset = parse_linemask_file(LMASK_DIR / "Ni" / "ni-lmask.txt")
    first = mset.masks[0]
    assert first.lam0 == pytest.approx(4811.9829, abs=1e-3)
    assert first.lmin == pytest.approx(4811.69, abs=1e-2)
    assert first.lmax == pytest.approx(4812.13, abs=1e-2)


def test_angstrom_file_unchanged():
    mset = parse_linemask_file(LMASK_DIR / "Fe" / "fe-lmask.txt")
    first = mset.masks[0]
    assert first.lam0 == pytest.approx(4808.148, abs=1e-3)
    assert first.lmin == pytest.approx(4807.97, abs=1e-2)
    assert first.lmax == pytest.approx(4808.33, abs=1e-2)


def test_two_column_segment_file():
    # Co/co-seg.txt is a headered 'wave_base wave_top' segment table: no centre, just bounds
    mset = parse_linemask_file(LMASK_DIR / "Co" / "co-seg.txt")
    assert mset.masks
    for m in mset.masks:
        assert math.isnan(m.lam0)
        assert m.lmin < m.lmax


def test_inline_comment_and_species_code():
    # Ba rows carry a species code (54.02) and one has a trailing '; blend with Fe 1' comment
    mset = parse_linemask_file(LMASK_DIR / "Ba" / "ba-lmask.txt")
    first = mset.masks[0]
    assert first.lam0 == pytest.approx(4934.076, abs=1e-3)
    assert first.lmin == pytest.approx(4933.945, abs=1e-3)
    assert first.lmax == pytest.approx(4934.500, abs=1e-3)
    assert "blend" not in first.note.lower()   # inline comment stripped, not folded into note


def test_multiple_masks_may_share_centre():
    # H-lmask.txt has three windows on the same Hβ centre (4861.9 Å)
    mset = parse_linemask_file(LMASK_DIR / "H" / "H-lmask.txt")
    centres = [m.lam0 for m in mset.masks]
    assert centres.count(pytest.approx(4861.9, abs=1e-3)) >= 2


def test_pad_controls_synth_window():
    mset = parse_linemask_file(LMASK_DIR / "Fe" / "fe-lmask.txt", pad=0.25)
    m = mset.masks[0]
    assert m.smin == pytest.approx(m.lmin - 0.25, abs=1e-9)
    assert m.smax == pytest.approx(m.lmax + 0.25, abs=1e-9)
    mset.repad(1.0)
    m = mset.masks[0]
    assert m.smin == pytest.approx(m.lmin - 1.0, abs=1e-9)


@pytest.mark.parametrize("path", ALL_FILES, ids=lambda p: str(p.relative_to(LMASK_DIR)))
def test_round_trip(tmp_path, path):
    mset = parse_linemask_file(path)
    if not mset.masks:
        pytest.skip("file has no active masks")
    out = tmp_path / "out.txt"
    write_linemask_file(out, mset)
    again = parse_linemask_file(out)
    assert len(again.masks) == len(mset.masks)
    for a, b in zip(mset.masks, again.masks):
        assert a.lmin == pytest.approx(b.lmin, abs=1e-3)
        assert a.lmax == pytest.approx(b.lmax, abs=1e-3)


def test_renorm_knots_are_colinear_and_span_window():
    m = LineMask(lam0=5000.0, lmin=4999.5, lmax=5000.5)
    m.pad_synth(0.5)
    m.renorm = True
    m.cont_pts = ((4999.5, 1.00), (5000.5, 1.10))   # a tilted local continuum
    knots = m.renorm_knots()
    assert knots is not None and len(knots) == 4
    xs = [k[0] for k in knots]
    assert xs[0] == pytest.approx(m.smin) and xs[-1] == pytest.approx(m.smax)
    # all four knots colinear ⇒ zero curvature
    xs = np.array([k[0] for k in knots])
    ys = np.array([k[1] for k in knots])
    slope = (ys[1] - ys[0]) / (xs[1] - xs[0])
    assert np.allclose(ys, ys[0] + slope * (xs - xs[0]), atol=1e-9)


class _StubCtx:
    """Minimal context: just the event bus MaskStore needs."""
    def __init__(self):
        self.events = []

    def emit(self, event, **data):
        self.events.append(event)


def test_maskstore_load_edit_and_broadcast():
    from maskstore import MaskStore
    ctx = _StubCtx()
    store = MaskStore(ctx)
    store.load_file(LMASK_DIR / "Fe" / "fe-lmask.txt")
    assert store.masks and store.active == 0
    assert "masks-changed" in ctx.events
    # resize the active fit window; synth window follows the pad
    store.set_pad(0.3)
    store.set_fit(0, 5000.0, 5000.6)
    m = store.active_mask()
    assert m.lmin == pytest.approx(5000.0) and m.lmax == pytest.approx(5000.6)
    assert m.smin == pytest.approx(5000.0 - 0.3) and m.smax == pytest.approx(5000.6 + 0.3)
    # add / remove keep the active index valid
    n = len(store.masks)
    idx = store.add(4000.0, 4000.5, note="test")
    assert idx == n and store.active == n
    store.remove(idx)
    assert len(store.masks) == n and store.active <= n - 1


def test_maskstore_snapshot_round_trip():
    from maskstore import MaskStore
    store = MaskStore(_StubCtx())
    store.load_file(LMASK_DIR / "Ni" / "ni-lmask.txt")
    store.set_renorm(0, True)
    snap = store.snapshot()
    store2 = MaskStore(_StubCtx())
    from tumult.linemask import LineMaskSet
    store2.restore([LineMaskSet.from_dict(d) for d in snap["sets"]],
                   snap["current"], snap["active"])
    assert len(store2.masks) == len(store.masks)
    assert store2.masks[0].renorm and store2.active == store.active


class _Raw:
    def __init__(self, x, y):
        self._x, self._y = np.asarray(x, float), np.asarray(y, float)

    @property
    def values(self):
        return self._x, self._y


def test_seed_cont_markers_fits_sloped_continuum_over_core():
    from linemask import seed_cont_markers
    x = np.linspace(5000.0, 5001.0, 201)
    cont = 1.0 + 0.5 * (x - 5000.0)          # a sloped continuum
    y = cont.copy()
    y[(x > 5000.45) & (x < 5000.55)] -= 0.4  # an absorption core the fit must ignore
    m = LineMask(lam0=5000.5, lmin=5000.4, lmax=5000.6)
    m.pad_synth(0.3)
    (x1, y1), (x2, y2) = seed_cont_markers(_Raw(x, y), m, pct=80)
    assert y1 == pytest.approx(1.0 + 0.5 * (x1 - 5000.0), abs=0.05)
    assert y2 == pytest.approx(1.0 + 0.5 * (x2 - 5000.0), abs=0.05)
    assert y2 > y1                           # positive continuum slope preserved


def test_linemask_serialisation_round_trip():
    m = LineMask(lam0=5000.0, lmin=4999.5, lmax=5000.5, note="Fe I", renorm=True,
                 cont_pts=((4999.5, 1.0), (5000.5, 1.1)))
    m.pad_synth(0.4)
    d = m.to_dict()
    m2 = LineMask.from_dict(d)
    assert m2.lam0 == m.lam0 and m2.note == m.note and m2.renorm
    assert m2.cont_pts == m.cont_pts
