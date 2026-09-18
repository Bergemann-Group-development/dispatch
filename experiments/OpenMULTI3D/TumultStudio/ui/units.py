"""Display units: the selectable axis units (label -> astropy unit), the conversions between
them, and the native unit the line list's wavelengths are converted from."""
import astropy.units as u

# native unit of the master line list's ``lam0`` column (:mod:`tumult.linelist`), which the
# line navigator and the line annotations convert to the canvas's display unit
LAM0_UNIT = u.AA

X_UNITS = {
    "Å": u.AA,
    "nm": u.nm,
    "µm": u.um,
    "eV": u.eV,
}

Y_UNITS = {
    "normalized": u.dimensionless_unscaled,
    "erg/s/cm²/Å": u.erg / (u.s * u.cm**2 * u.AA),
    "erg/s/cm²/nm": u.erg / (u.s * u.cm**2 * u.nm),
}


def label_for(unit, mapping, default):
    """Reverse-lookup the menu label for an astropy unit."""
    for label, candidate in mapping.items():
        if candidate == unit:
            return label
    return default


def to_value(value, from_unit, to_unit):
    """Convert a value (scalar or array) between spectral units.

    Uses the ``spectral()`` equivalency so wavelength ↔ frequency ↔ energy all work.
    Returns a plain number/array (the ``.value``), not a Quantity.
    """
    return (value * from_unit).to(to_unit, equivalencies=u.spectral()).value


def to_range(lo, hi, from_unit, to_unit):
    """Convert an ``(lo, hi)`` interval, returned ascending.

    A unit flip (e.g. wavelength → energy) reverses the order, so the result is always
    sorted low-to-high.
    """
    a = to_value(lo, from_unit, to_unit)
    b = to_value(hi, from_unit, to_unit)
    return (a, b) if a <= b else (b, a)
