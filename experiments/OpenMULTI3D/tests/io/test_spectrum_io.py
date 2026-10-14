"""Pin the observation readers moved out of the Bokeh app into :mod:`tumult.spectrum_io`.

Observation files arrive in whatever shape their survey wrote them, so the delimiter and
header size are sniffed rather than configured. These fix that sniffing against the shapes
seen in practice — bare columns, comment banners, a leading point count, a title line, and
the comma/tab separated variants — plus the ascending-wavelength guarantee ``Spectrum.crop``
relies on for its ``searchsorted``.
"""
import numpy as np
import pytest

from tumult import spectrum_io as sio

LAM = np.linspace(5000.0, 5001.0, 11)
FLUX = 1.0 - 0.3 * np.exp(-((LAM - 5000.5) / 0.1) ** 2)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def _rows(sep=" ", lam=LAM, flux=FLUX):
    return "\n".join(f"{a}{sep}{b}" for a, b in zip(lam, flux))


@pytest.mark.parametrize("name,text", [
    ("plain.txt", _rows()),
    ("comment.txt", "# wavelength flux\n# from somewhere\n" + _rows()),
    ("count.txt", f"{len(LAM)}\n" + _rows(sep="  ")),
    ("title.txt", "SOME SURVEY SPECTRUM\n" + _rows()),
    ("csv.csv", _rows(sep=",")),
    ("tab.txt", _rows(sep="\t")),
    ("semi.txt", _rows(sep=";")),
])
def test_sniffs_delimiter_and_header(tmp_path, name, text):
    spec = sio.load_spectrum(_write(tmp_path, name, text))
    x, y = spec.values
    assert len(x) == len(LAM)
    assert np.allclose(x, LAM)
    assert np.allclose(y, FLUX)


def test_descending_input_is_flipped_ascending(tmp_path):
    """crop's searchsorted needs ascending x, whatever order the file was written in."""
    text = _rows(lam=LAM[::-1], flux=FLUX[::-1])
    x, y = sio.load_spectrum(_write(tmp_path, "desc.txt", text)).values
    assert x[0] < x[-1]
    assert np.allclose(x, LAM)
    assert np.allclose(y, FLUX)


def test_unparseable_rows_are_dropped_not_fatal(tmp_path):
    """One bad row must not abort the load of an otherwise good file."""
    rows = _rows().splitlines()
    rows.insert(5, "NaN  ---")
    x, _y = sio.load_spectrum(_write(tmp_path, "junk.txt", "\n".join(rows))).values
    assert len(x) == len(LAM)


def test_blank_file_yields_no_points(tmp_path):
    with pytest.raises(ValueError):
        sio.load_spectrum(_write(tmp_path, "empty.txt", "# only a comment\n"))


def test_is_run_needs_params_nml(tmp_path):
    assert not sio.is_run(tmp_path)
    (tmp_path / "params.nml").write_text("&io_params /\n")
    assert sio.is_run(tmp_path)


def test_load_file_dispatches_on_suffix(tmp_path):
    """A .txt path goes through the text reader; the FITS branch is suffix-selected."""
    p = _write(tmp_path, "obs.txt", _rows())
    spec = sio.load_file(p)
    assert np.allclose(spec.values[0], LAM)
