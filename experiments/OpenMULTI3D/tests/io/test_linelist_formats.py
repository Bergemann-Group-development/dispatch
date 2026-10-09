"""Line lists in the tumult and Turbospectrum layouts read into the same frame.

The layout is told per block (tumult.linelist._blocks), and the isotope label is
Turbospectrum's mass part in both (tumult.isotopes).
"""
import numpy as np
import pandas as pd

from tumult.isotopes import iso_token, mass_label, ts_species, ts_token
from tumult.linelist import LineList, override_lines

TS = """'  26.056   '    1    2
'Fe I    NLTE'
  5161.738 4.294 -5.246 -7.510 7.0 2.09E+08 0.000 's' 'p' 0.0 1.0 'b1D t5P*'  102 422  'a' 'b'
  5162.272 4.178 0.020 798.279 11.0 2.40E+08 'p' 'd' 0.0 1.0  91 378  'none' 'none'
'0607.013014 '    1    1
'CN'
  5163.000 0.263 -5.828 .00 6.5 1.000E+00 'X' 'X'   0.0    1.0 'x'
"""

TUMULT = """'Fe I'
'56' 2
  5161.738 4.294 -5.246 -7.510 2.09E+08 0.000 102 422
  5162.272 4.178 0.020 798.279 2.40E+08 0 0 0
'CN'
'013014' 1
  5163.000 0.263 -5.828 0 1.000E+00 0 0 0
"""

COLUMNS = ['lam0', 'evi', 'loggf', 'gvdw', 'grad', 'gstk', 'i', 'j']


def test_both_layouts_read_into_the_same_lines(tmp_path):
    (tmp_path / 'ts').write_text(TS)
    (tmp_path / 'tumult').write_text(TUMULT)
    ts = LineList.from_ts_format(tmp_path / 'ts')
    tumult = LineList.from_tumult_format(tmp_path / 'tumult')
    assert list(ts.index.droplevel(-1).unique()) == [('CN', 13014), ('Fe I', 56)]
    # the second Fe I record has no comment, so its indices are not read (as in
    # Turbospectrum, which reads such a block without them)
    pd.testing.assert_frame_equal(ts[COLUMNS], tumult[COLUMNS], check_dtype=False)
    assert list(ts.levlo) == ['X', 's', 'p']


def test_a_file_may_mix_the_layouts(tmp_path):
    (tmp_path / 'mixed').write_text(TUMULT + TS)
    mixed = LineList.from_tumult_format(tmp_path / 'mixed')
    assert len(mixed) == 6
    assert set(mixed.index.droplevel(-1)) == {('CN', 13014), ('Fe I', 56)}


def test_isotope_labels_follow_turbospectrum():
    assert ts_species('0106.001012', 1) == ('CH', 12)          # H always natural
    assert ts_species('0606.013012', 1) == ('C2', 12013)       # equal atoms ascending
    assert ts_species('26.000', 2) == ('Fe II', 0)
    assert mass_label([(8, 16), (1, 1), (1, 1)]) == 16
    assert iso_token('CH', 12) == '000012'
    assert iso_token('Fe I', 56) == '56'
    assert ts_token('CN', 12014) == '0607.012014'
    assert ts_token('Fe I', 0) == '26.000'


def test_tumult_writer_keeps_the_molecular_label(tmp_path):
    (tmp_path / 'ts').write_text(TS)
    ts = LineList.from_ts_format(tmp_path / 'ts')
    ts.reset_index(level=[0, 1]).pipe(LineList).write_to_tumult_format(tmp_path / 'out')
    text = (tmp_path / 'out').read_text()
    assert "'013014' 1" in text
    back = LineList.from_tumult_format(tmp_path / 'out')
    np.testing.assert_allclose(back.lam0, ts.lam0)
    assert list(back.index.droplevel(-1).unique()) == [('CN', 13014), ('Fe I', 56)]


def test_override_edits_only_the_matched_record(tmp_path):
    for name, text, line in (('ts', TS, dict(ion='Fe I', lam0=5162.28, loggf=-0.5)),
                             ('tumult', TUMULT, dict(ion='Fe I', lam0=5162.28, loggf=-0.5))):
        (tmp_path / name).write_text(text)
        out = override_lines([line], tmp_path / f'{name}.out', base=tmp_path / name)
        old, new = text.splitlines(), out.read_text().splitlines()
        changed = [k for k, (o, n) in enumerate(zip(old, new)) if o != n]
        assert len(old) == len(new) and len(changed) == 1
        fe = LineList.from_tumult_format(out).loc[('Fe I', 56)]
        # the record keeps its wavelength and the columns the line leaves out
        assert list(fe.loggf.round(3)) == [-5.246, -0.5]
        assert list(fe.lam0) == [5161.738, 5162.272] and fe.gvdw.iloc[1] == 798.279
        if name == 'ts':
            # the Turbospectrum record keeps its level designations and comment
            k = changed[0]
            assert new[k].split("'", 1)[1] == old[k].split("'", 1)[1]


def test_override_rejects_an_ambiguous_record(tmp_path):
    (tmp_path / 'base').write_text(TS + TUMULT)
    try:
        override_lines([dict(ion='CN', lam0=5163.0, gvdw=250.0)], tmp_path / 'out',
                       base=tmp_path / 'base')
    except ValueError as err:
        assert 'equally near' in str(err)
    else:
        raise AssertionError("a line with two equally near records was accepted")


def test_override_adds_a_line_without_record(tmp_path):
    (tmp_path / 'base').write_text(TUMULT)
    out = override_lines([dict(ion='Fe I', iso=56, lam0=5170.0, evi=1.0, loggf=-2.0)],
                         tmp_path / 'out', base=tmp_path / 'base')
    fe = LineList.from_tumult_format(out).loc[('Fe I', 56)]
    assert list(fe.lam0) == [5161.738, 5162.272, 5170.0]
    try:
        override_lines([dict(ion='Fe I', lam0=5170.0, loggf=-2.0)], tmp_path / 'bad',
                       base=tmp_path / 'base')
    except ValueError as err:
        assert "needs 'evi'" in str(err)
    else:
        raise AssertionError("a new line without evi was accepted")
