"""Parser and round-trip tests for the LineMask model (``tumult/linemask.py``).

Every file under ``input/linemasks/`` must parse into sane, ascending-bound masks with
wavelengths normalised to Å, and the canonical writer must round-trip them. Those files
are local input data that the repository does not carry.
"""
import math
from pathlib import Path

import numpy as np
import pytest

from tumult.linemask import LineMaskSet, parse_linemask_file, write_linemask_file

pytestmark = pytest.mark.requires_input_data

ROOT = Path(__file__).resolve().parents[2]
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
