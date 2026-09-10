"""Reading observed spectra off disk: two-column text/CSV, FITS, and run folders.

The data model is :class:`tumult.Spectrum` (crop, slicing, units, arithmetic, get_weq,
fit, ...), reused as-is; what this module adds is the loading it does not do. Text
delimiters and header sizes are sniffed rather than configured, because observation files
arrive in whatever shape their survey wrote them.

Kept apart from :mod:`tumult.spectrum` so the data model stays free of pandas and
astropy.io.fits. Windowed decimation is just ``spectrum.crop(xmin=..., xmax=...)[::step]``.
"""
from __future__ import annotations

from pathlib import Path

import astropy.io.fits as fits
import astropy.units as u
import numpy as np
import pandas as pd

from .spectrum import Spectrum  # shared data model — re-exported, never reimplemented
from .run import Run as _Run  # run reader (Run.read / get_spectrum)

__all__ = ["Spectrum", "load_spectrum", "load_file", "is_run", "load_run"]


def _sniff(path):
    """Sniff ``(sep, skiprows)`` for a two-column spectrum file.

    Scans from the top and infers the header size: every leading line that is
    blank, a ``#`` comment, or otherwise not a data row is counted toward the
    header. A line only counts as the data start when it splits into **at least
    two float-parseable fields** — so a lone leading integer (a point count) is
    treated as header, not data. The first such data line fixes the delimiter
    (the first single-char ``, \\t ; |`` present, else any whitespace run) and
    its index becomes ``skiprows``.
    """
    seps = (",", "\t", ";", "|")
    skiprows = 0
    sep = r"\s+"
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                skiprows += 1
                continue
            this_sep = next((c for c in seps if c in s), r"\s+")
            fields = s.split(this_sep) if this_sep != r"\s+" else s.split()
            floats = 0
            for tok in fields[:2]:
                try:
                    float(tok)
                except ValueError:
                    break
                floats += 1
            if floats >= 2:
                return this_sep, skiprows
            skiprows += 1
    return sep, skiprows


def load_spectrum(path, **read_kwargs) -> Spectrum:
    """Load a two-column (wavelength, flux) text/CSV spectrum.

    The delimiter (whitespace or comma/tab/semicolon) and the header size (any
    number of leading comment/title/count lines) are auto-detected; pass explicit
    ``pd.read_csv`` kwargs to override.
    """
    path = Path(path)
    sep, skiprows = _sniff(path)
    kwargs = dict(sep=sep, header=None, skiprows=skiprows, comment="#",
                  engine="python" if sep == r"\s+" else "c")
    kwargs.update(read_kwargs)
    df = pd.read_csv(path, **kwargs)
    # coerce non-numeric junk (stray count lines, embedded NULs, ...) to NaN and
    # drop those rows rather than letting a single bad row abort the whole load
    x = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy()
    y = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy()
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if x.size and x[0] > x[-1]:  # keep ascending so crop's searchsorted is happy
        x, y = x[::-1], y[::-1]
    # no units in a plain text file -> wavelength assumed Å, flux assumed normalized
    return Spectrum(x, y, yunit=u.dimensionless_unscaled)


def _cards_frame(header):
    """A ``(key, value, comment)`` DataFrame of a FITS header's cards."""
    return pd.DataFrame(data=np.array(header.cards),
                        columns=["key", "value", "comment"])


def _read_4most_fits(file):
    """Read a 4MOST/QMOST FITS spectrum: WAVE/FLUX columns in HDU[1], units from
    its TUNIT keywords; header cards from HDUs 0 and 1."""
    with fits.open(file) as hdul:
        header0, header1 = hdul[0].header, hdul[1].header
        data = hdul[1].data
        meta = pd.concat([_cards_frame(header0), _cards_frame(header1)])
        spec = Spectrum(data.WAVE.squeeze(), data.FLUX.squeeze(),
                        xunit=header1["TUNIT1"], yunit=header1["TUNIT2"])
    return spec, meta


def _tbdfits(file, ext=0, survey=None):
    """Read a generic WCS-sampled FITS spectrum: flux in the data, wavelength
    built from ``CRVAL1`` + ``CD1_1``/``CDELT1`` (log-sampled for survey ``apo``)."""
    with fits.open(file) as hdul:
        header = hdul[ext].header
        flux = hdul[ext].data

    cr = header.get("CRVAL1")
    if cr is None:
        raise ValueError("CRVAL1 not found in the header.")
    cdelt = header.get("CD1_1") or header.get("CDELT1")
    if cdelt is None:
        raise ValueError("CDELT1/CD1_1 not found in the header.")

    lam = np.arange(len(flux)) * cdelt + cr
    if survey == "apo":
        lam = 10 ** lam
    return Spectrum(lam, flux), _cards_frame(header)


def _fits_meta(df):
    """A JSON-serialisable Info-tab payload from a FITS header-card DataFrame."""
    cards = [[str(k), str(v), str(c)]
             for k, v, c in df[["key", "value", "comment"]].itertuples(index=False)
             if str(k).strip()]
    return {"type": "fits", "cards": cards}


def load_file(path) -> Spectrum:
    """Load a spectrum from a text or FITS file.

    For FITS, the header cards are attached to the returned spectrum as ``.meta`` (a
    structured dict the Info tab renders); text files carry no metadata.
    """
    path = Path(path)
    if path.suffix.lower() in (".fits", ".fit"):
        try:
            spec, meta = _read_4most_fits(str(path))
        except Exception:
            spec, meta = _tbdfits(str(path))
        spec.meta = _fits_meta(meta)
        return spec
    return load_spectrum(path)


def is_run(folder) -> bool:
    """A tumult run folder is identified by a ``params.nml`` directly inside it."""
    return (Path(folder) / "params.nml").is_file()


def load_run(folder) -> Spectrum:
    """Load a tumult run's continuum-normalized spectrum (flux / continuum).

    The result is dimensionless (≈0–1), so it overlays directly on the normalized
    objects the canvas works with.
    """
    run = _Run.read(Path(folder))
    return run.get_spectrum() / run.get_spectrum(cont=True)
