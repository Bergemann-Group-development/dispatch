"""Spectrum loading for TumultStudio2.

The data model is tumult's `Spectrum` (crop, slicing, units, arithmetic,
get_weq, fit, ...), reused as-is. The only thing it lacks for the viewer is
loading from a text file; windowed decimation is just
``spectrum.crop(xmin=..., xmax=...)[::step]``.
"""
from __future__ import annotations

from pathlib import Path

import astropy.units as u
import pandas as pd

from tumult import Spectrum  # shared data model — re-exported, never reimplemented
from tumult.tumult import tumult as _Run  # run reader (tumult.read / get_spectrum)

__all__ = ["Spectrum", "load_spectrum", "load_file", "is_run", "load_run"]


def _sniff(path):
    """Sniff ``(sep, header)`` for a two-column spectrum file.

    pandas' own ``sep=None`` sniffer (csv.Sniffer) handles comma/tab/etc. but mangles
    whitespace-delimited files, so: use the first single-char delimiter present, else
    fall back to any whitespace run; treat the first line as a header iff its first
    field isn't numeric.
    """
    with open(path) as f:
        first = next((s for s in (ln.strip() for ln in f)
                      if s and not s.startswith("#")), "")
    sep = next((c for c in (",", "\t", ";", "|") if c in first), r"\s+")
    tokens = first.split(sep) if sep != r"\s+" else first.split()
    try:
        float(tokens[0])
        header = None
    except (ValueError, IndexError):
        header = 0
    return sep, header


def load_spectrum(path, **read_kwargs) -> Spectrum:
    """Load a two-column (wavelength, flux) text/CSV spectrum.

    The delimiter (whitespace or comma/tab/semicolon) and whether there's a header
    row are auto-detected; pass explicit ``pd.read_csv`` kwargs to override.
    """
    path = Path(path)
    sep, header = _sniff(path)
    kwargs = dict(sep=sep, header=header, comment="#",
                  engine="python" if sep == r"\s+" else "c")
    kwargs.update(read_kwargs)
    df = pd.read_csv(path, **kwargs)
    x = df.iloc[:, 0].to_numpy(dtype=float)
    y = df.iloc[:, 1].to_numpy(dtype=float)
    if x.size and x[0] > x[-1]:  # keep ascending so crop's searchsorted is happy
        x, y = x[::-1], y[::-1]
    # no units in a plain text file -> wavelength assumed Å, flux assumed normalized
    return Spectrum(x, y, yunit=u.dimensionless_unscaled)


def _fits_meta(df):
    """A JSON-serialisable Info-tab payload from a FITS header-card DataFrame."""
    cards = [[str(k), str(v), str(c)]
             for k, v, c in df[["key", "value", "comment"]].itertuples(index=False)
             if str(k).strip()]
    return {"type": "fits", "cards": cards}


def load_file(path) -> Spectrum:
    """Load a spectrum from a text or FITS file (FITS via the survey readers).

    For FITS, the header cards are attached to the returned spectrum as ``.meta`` (a
    structured dict the Info tab renders); text files carry no metadata.
    """
    path = Path(path)
    if path.suffix.lower() in (".fits", ".fit"):
        # Bokeh restores sys.path after app startup, dropping this app directory,
        # so re-add it before importing sibling modules from a callback.
        import sys
        here = str(Path(__file__).parent)
        if here not in sys.path:
            sys.path.insert(0, here)
        import surveys
        try:
            spec, meta = surveys.read_4most_fits(str(path))
        except Exception:
            spec, meta = surveys.tbdfits(str(path))
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
