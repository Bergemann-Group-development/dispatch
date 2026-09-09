"""The handful of atmospheres a strength ranking is ever actually computed on.

Why a table and not the grid
----------------------------
A strength is a maximum over depth of line-centre over continuum extinction, and
:func:`linestrength.compute_band` needs exactly two columns of an atmosphere to get it —
T and rho, at eight depths spread over the model's whole optical-depth range (the sampling
argument is in :func:`linestrength.read_marcs`'s ancestor; eight reproduces all 56 layers to
Spearman 1.0000). Eight pairs of numbers per star is a table, not a file format: the MARCS
grid it came from is fifteen thousand models on disk that have to be found, name-parsed and
opened, all to recover sixteen floats.

So the sixteen floats are here, extracted once from the MARCS models named in each entry, and
the grid is not a dependency any more. What is lost is the ability to land on the *exact*
model of an arbitrary star. That was never worth much: the ranking is what the strengths are
for, and it moves slowly with the parameters — which is the same reason the old code was
content with the nearest grid point rather than an interpolation.

Choosing between them
---------------------
:func:`nearest` picks by the same scaled distance the grid lookup used — 1000 K of Teff,
1 dex of gravity and 1 dex of metallicity weighted alike — so a spectrum that states its own
parameters still lands on a sensible atmosphere without anyone choosing one. The names are
there for the spectra that state nothing, which is most observed ones: a class is something
you can pick off a dropdown from what you know about the star, where a Teff in Kelvin is a
number you have to invent.

The table is deliberately coarse. Nine classes spanning 2500-6500 K and log g 0-5 cover the
cool half of the HR diagram where a line list has anything to say; adding entries is a matter
of naming a model and running the extractor in this module's docstring:

    T, rho = linestrength.read_marcs(path)   # with NDEPTH layers
"""
from __future__ import annotations

import numpy as np


class RefStar:
    """One reference atmosphere: what to call it, and the T/rho HOP needs.

    ``teff``/``logg``/``feh`` are the model's own parameters. They are what :func:`nearest`
    matches on and what the panel footer shows; nothing recomputes an atmosphere from them.
    """

    __slots__ = ("key", "label", "teff", "logg", "feh", "model", "T", "rho")

    def __init__(self, key, label, teff, logg, feh, model, T, rho):
        self.key = key
        self.label = label
        self.teff = float(teff)
        self.logg = float(logg)
        self.feh = float(feh)
        self.model = model          # the MARCS file the two arrays were extracted from
        self.T = np.asarray(T, float)
        self.rho = np.asarray(rho, float)

    @property
    def summary(self):
        """``Teff 5750, log g 4.50, [Fe/H] +0.00`` — the parameters, for a footer."""
        return f"Teff {self.teff:.0f}, log g {self.logg:.2f}, [Fe/H] {self.feh:+.2f}"

    def __repr__(self):
        return f"<RefStar {self.key} {self.summary}>"


# Extracted from the named MARCS models at NDEPTH=8 depths, evenly spread over the model's
# own depth index — see the module docstring. T in K, rho in g/cm3, surface first.
REFERENCE_STARS = [
    RefStar("gdwarf", "G dwarf (solar)", 5750, 4.5, +0.00,
            "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00",
            [4045.9, 4395.3, 4616.0, 4815.4, 5124.9, 6017.7, 8057.7, 9909.1],
            [1.085e-09, 7.663e-09, 2.550e-08, 6.602e-08,
             1.466e-07, 3.037e-07, 3.251e-07, 3.625e-07]),
    RefStar("fdwarf", "F dwarf", 6500, 4.0, +0.00,
            "p6500_g+4.0_m0.0_t01_st_z+0.00_a+0.00",
            [4606.5, 4968.3, 5218.4, 5442.5, 5768.4, 6766.6, 8996.9, 10749.0],
            [4.608e-10, 3.318e-09, 1.047e-08, 2.584e-08,
             5.343e-08, 8.696e-08, 7.784e-08, 7.483e-08]),
    RefStar("kdwarf", "K dwarf", 4500, 4.5, +0.00,
            "p4500_g+4.5_m0.0_t01_st_z+0.00_a+0.00",
            [3116.7, 3470.0, 3644.5, 3802.1, 4068.8, 4796.0, 6505.1, 8634.2],
            [3.724e-09, 2.255e-08, 6.972e-08, 1.740e-07,
             3.517e-07, 5.964e-07, 8.133e-07, 1.031e-06]),
    RefStar("mdwarf", "M dwarf", 3500, 5.0, +0.00,
            "p3500_g+5.0_m0.0_t01_st_z+0.00_a+0.00",
            [2607.6, 2709.5, 2831.7, 3063.3, 3358.1, 3742.5, 4155.9, 4932.7],
            [2.572e-09, 2.762e-08, 1.535e-07, 5.438e-07,
             1.394e-06, 3.571e-06, 7.678e-06, 2.132e-05]),
    RefStar("turnoff", "Turn-off (metal-poor)", 6250, 4.0, -2.00,
            "p6250_g+4.0_m0.0_t01_st_z-2.00_a+0.40",
            [4638.6, 4925.3, 5046.9, 5160.8, 5440.5, 6443.4, 8414.4, 10047.0],
            [6.331e-10, 5.314e-09, 1.918e-08, 5.508e-08,
             1.190e-07, 1.782e-07, 1.592e-07, 1.621e-07]),
    RefStar("clump", "Red clump", 4750, 2.5, +0.00,
            "s4750_g+2.5_m1.0_t01_st_z+0.00_a+0.00",
            [3177.7, 3570.2, 3797.7, 3996.6, 4280.9, 5034.0, 6900.1, 8970.9],
            [1.432e-10, 1.028e-09, 3.112e-09, 7.611e-09,
             1.605e-08, 3.439e-08, 4.394e-08, 3.887e-08]),
    RefStar("rgb", "Red giant branch", 4500, 2.0, +0.00,
            "s4500_g+2.0_m1.0_t01_st_z+0.00_a+0.00",
            [2954.6, 3364.5, 3590.6, 3785.5, 4071.5, 4797.4, 6538.1, 8748.8],
            [1.065e-10, 7.370e-10, 2.112e-09, 4.961e-09,
             9.941e-09, 1.996e-08, 2.735e-08, 2.307e-08]),
    RefStar("rgbmp", "Metal-poor giant", 4500, 1.5, -2.00,
            "s4500_g+1.5_m1.0_t02_st_z-2.00_a+0.40",
            [3001.5, 3482.1, 3656.9, 3803.6, 4072.5, 4788.6, 6321.7, 7823.8],
            [1.795e-10, 1.581e-09, 5.485e-09, 1.439e-08,
             3.058e-08, 4.959e-08, 4.401e-08, 3.944e-08]),
    RefStar("rsg", "Red supergiant", 3800, 0.0, +0.00,
            "s3800_g+0.0_m15._t05_st_z+0.00_a+0.00",
            [2464.7, 2794.9, 3026.0, 3231.6, 3508.1, 4155.1, 5580.3, 8145.1],
            [4.004e-12, 5.879e-11, 2.393e-10, 6.074e-10,
             1.161e-09, 1.927e-09, 3.051e-09, 2.232e-09]),
]

BY_KEY = {s.key: s for s in REFERENCE_STARS}
DEFAULT_KEY = "gdwarf"
AUTO_KEY = "auto"                # follow the spectrum's own parameters, when it has any


def options(auto_label="Auto (from spectrum)"):
    """``[(key, label)]`` for a select field, ``auto`` first."""
    return [(AUTO_KEY, auto_label)] + [(s.key, s.label) for s in REFERENCE_STARS]


def get(key):
    """The star named ``key``, or the default one. ``auto`` is not a star — it resolves
    against a spectrum's parameters, so callers wanting one go through
    :meth:`linestrength.LineStrengths.star_for`."""
    return BY_KEY.get(str(key), BY_KEY[DEFAULT_KEY])


def nearest(teff, logg, feh):
    """The table entry closest to ``(teff, logg, feh)``.

    Distances are scaled so the three are comparable: 1000 K of Teff, 1 dex of gravity and
    1 dex of metallicity count the same — the weighting the MARCS grid lookup used. Teff
    dominates in practice, which is what lets a nine-entry table stand in for a grid.
    """
    t = np.array([s.teff for s in REFERENCE_STARS])
    g = np.array([s.logg for s in REFERENCE_STARS])
    z = np.array([s.feh for s in REFERENCE_STARS])
    d = ((t - float(teff)) / 1000.0) ** 2 + (g - float(logg)) ** 2 + (z - float(feh)) ** 2
    return REFERENCE_STARS[int(np.argmin(d))]
