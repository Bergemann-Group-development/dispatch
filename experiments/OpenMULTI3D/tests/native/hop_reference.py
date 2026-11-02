"""The HOP EOS and opacity states the native references are computed on.

Inputs are the tracked defaults of ``hopacity.HopEOS`` (``input/absdat/TS_absdat.txt``,
``input/abund/abund_magg22.txt``, ``input/isotopes/isotopes_blp.txt``,
``input/molec/molec_data_BC.txt``) and ``input/linelists/master_line.txt``, through the
``libhop.so`` that ``make -C HOP/standalone`` copies into ``HOP/hopacity/``.
"""
from pathlib import Path

import numpy as np

PROJECT_DIR = Path(__file__).resolve().parents[2]
HOP_LIB = PROJECT_DIR / "HOP" / "hopacity" / "libhop.so"
HOP_LINELIST = PROJECT_DIR / "input" / "linelists" / "master_line.txt"

# a coarse T x rho node grid (4 x 6 = 24 cells) and a narrow 10 AA window: small enough to
# solve in seconds, wide enough to cross the Mg b lines and a good deal of Fe I
HOP_DLAM = 10.0        # line-wing search radius [AA]
HOP_VMIC = 1.0         # [km/s]; must be pinned, it moves the line opacity a lot


def hop_grid():
    """The (temp [K], rho [g/cm^3]) node grid, both shaped (ntemp, nrho)."""
    from hopacity import logstep

    tab_temp = logstep(3500.0, 12000.0, 0.35)
    tab_rho = logstep(1e-10, 1e-7, 1.2)
    return np.meshgrid(tab_temp, tab_rho, indexing="ij")


def hop_wavelengths():
    """The vacuum wavelength grid [AA] of the opacity reference."""
    from hopacity import wavelength_grid

    return wavelength_grid(5160.0, 5170.0, R_sampling=20000)


def hop_eos_reference():
    """Thermodynamic state, two species densities and the (T, ne) -> rho inversion."""
    from hopacity import HopEOS

    temp, rho = hop_grid()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho, species=["Fe I", "CH"])
        out = {
            "temp": hop.temp, "rho": hop.rho, "ne": hop.ne, "pe": hop.pe, "pg": hop.pg,
            "energy": hop.energy, "eh": hop.eh, "chi500": hop.chi500,
            "grph": np.array(hop.grph),
            "n_FeI": hop.species_density("Fe I"),
            "n_CH": hop.species_density("CH"),
        }
        # the (T, ne) branch of load_state must invert the (T, rho) branch
        hop.load_state(temp, ne=out["ne"])
        out["rho_from_ne"] = hop.rho
    return out


def hop_opacity_reference():
    """Continuum and line-list opacity [1/cm] and scattering on the node grid."""
    from hopacity import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        chi_c, scat_c = hop.opacity(lam)
        chi_l, scat_l = hop.opacity(lam, linelist=str(HOP_LINELIST), dlam=HOP_DLAM,
                                    vmic=HOP_VMIC, line_scattering=True)
        return {
            "lam": lam,
            "chi_continuum": chi_c, "scat_continuum": scat_c,
            "chi_linelist": chi_l, "scat_linelist": scat_l,
        }
