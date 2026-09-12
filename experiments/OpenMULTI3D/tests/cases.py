"""Deterministic characterization cases for the physics-critical numerics.

Each ``*_cases()`` function computes a flat ``{name: value}`` dict from the current code.
``gen_goldens.py`` dumps these to ``tests/golden/*.npz``; the test modules recompute them
and compare. Inputs are fixed (no randomness at measure time), so recomputation on the
same machine/env is expected to match to machine precision.
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# --- shared synthetic line -----------------------------------------------------------------
# Gaussian absorption line: amp=0.42, sig=0.11 AA at 5000.05 AA on a +-2 AA window.
# Analytic EW = amp*sqrt(2*pi)*sig = 0.11576... AA = 115.76 maa (window truncation negligible).
LINE_X0 = 5000.0
LINE_CENTER = 5000.05
LINE_AMP = 0.42
LINE_SIG = 0.11


def gauss_line(n=801, span=2.0):
    x = np.linspace(LINE_X0 - span, LINE_X0 + span, n)
    y = 1 - LINE_AMP * np.exp(-0.5 * ((x - LINE_CENTER) / LINE_SIG) ** 2)
    return x, y


# --- equivalent widths / profiles ----------------------------------------------------------

def ew_profile_cases():
    from tumult.profile import profile_names
    from tumult.utils import get_weq, fit_line_profile

    x, y = gauss_line()
    sigma = np.full(x.size, 1.0 / 80.0)  # SNR-80 per-pixel error model (no added noise)
    out = {}
    for name in sorted(profile_names()):
        prof = get_weq(x, y, x0=LINE_X0, profile=name, return_profile=True,
                       sigma=sigma, syst=0.02)
        # These dict keys are the stored golden keys and keep their original
        # spelling: renaming one silently misses its array in the .npz. The
        # attributes they read follow the lowercase unit-suffix convention.
        out[f"{name}_mAA"] = prof.maa
        out[f"{name}_err_phot_mAA"] = prof.err_phot_maa
        out[f"{name}_err_cont_mAA"] = prof.err_cont_maa
        out[f"{name}_rew"] = prof.rew
        out[f"{name}_rew_err"] = prof.rew_err
        if prof.fitted:
            out[f"{name}_fwhm"] = prof.fwhm
            out[f"{name}_params"] = np.asarray(prof.params, dtype=float)

    # the three public return conventions of get_weq (reduced EW default, width Quantity,
    # width + error tuple) must keep giving the same numbers through the refactor
    out["default_rew"] = get_weq(x, y, x0=LINE_X0)
    w, e = get_weq(x, y, x0=LINE_X0, target_unit=u.mAA, return_err=True, sigma=sigma)
    out["width_mAA"] = w.value
    out["width_err_mAA"] = e.value

    # direct fit entry point (pout, pcov) for the gaussian
    pout, pcov, _func = fit_line_profile(x, y, profile="gauss", x0=LINE_X0)
    out["fit_pout"] = np.asarray(pout, dtype=float)
    out["fit_pcov"] = np.asarray(pcov, dtype=float)
    return out


# --- broadening / convolution ---------------------------------------------------------------

def convolution_cases():
    from tumult.convolutions import combined_kernel, convert_res_to_vdop
    from tumult.spectrum import Spectrum

    out = {}
    kms = u.km / u.s
    out["res_to_vdop_50k_kms"] = convert_res_to_vdop(50000, None).to(kms).value
    out["res_to_vdop_50k_plus2_kms"] = convert_res_to_vdop(50000, 2.0 * kms).to(kms).value
    out["res_to_vdop_pair_kms"] = convert_res_to_vdop((80000, 20000), None).to(kms).value

    rd = 2e-6  # log-lambda step, dimensionless (v/c)
    px, py = combined_kernel(rd=rd, vdop=3.0 * kms, vrot=2.0 * kms, zeta_rt=1.5 * kms)
    out["kernel_px"] = np.asarray(px, dtype=float)
    out["kernel_py"] = np.asarray(py, dtype=float)

    # full crop + convolution round trip through Spectrum (plain-float kwargs, the GUI path)
    x, y = gauss_line()
    sp = Spectrum(x, y)
    conv = sp.crop(vdop=3.0, vrot=2.0, zeta_rt=1.5, vshift=0.5)
    out["crop_conv_x"] = np.asarray(conv.xdata.value, dtype=float)
    out["crop_conv_y"] = np.asarray(conv.ydata, dtype=float)
    res = sp.crop(res=50000)
    out["crop_res_y"] = np.asarray(res.ydata, dtype=float)
    return out


# --- air/vacuum wavelength conversion --------------------------------------------------------

AIRVAC_LINES = np.array([3933.66, 4554.033, 5167.32, 6562.79, 8542.09, 15000.0])  # AA (air-ish)


def airvac_cases():
    from tumult.auxfuncs import air2vac, vac2air

    grid = np.linspace(3000.0, 25000.0, 221)
    return {
        "vac2air_grid": vac2air(grid),
        "air2vac_grid": air2vac(grid),
        "vac2air_lines": vac2air(AIRVAC_LINES),
        "air2vac_lines": air2vac(AIRVAC_LINES),
    }


# --- abundance numerics (library + the GUI functions that move in Phase 2) -------------------

def _mk_line(a, sigma_a, phot, include=True):
    # continuum share derived exactly so sigma_a**2 = phot**2 + cont**2 holds to the bit
    cont = float(np.sqrt(sigma_a**2 - phot**2))
    return SimpleNamespace(a=a, sigma_a=sigma_a, sigma_a_phot=phot, sigma_a_cont=cont,
                           include=include)


AGG_LINES = [
    _mk_line(7.52, 0.05, 0.04),
    _mk_line(7.48, 0.08, 0.05),
    _mk_line(7.61, 0.04, 0.02),
    _mk_line(7.95, 0.30, 0.25),   # outlier, heavily down-weighted
    _mk_line(7.55, 0.06, 0.06, include=False),  # excluded by the user
]


def abundance_cases():
    # the GUI's `abundance` is a re-export of this module; the goldens pin the one
    # implementation both sides call.
    from tumult import abundance as gui_ab
    from tumult.abundance import aggregate, cog_invert, reduced_ew
    from tumult.spectrum import Spectrum

    out = {}
    # curve-of-growth inversion on a CoG that is exactly linear in (REW, A), so the REW-space
    # inversion recovers log10(50/55)/0.62
    a_grid = np.linspace(-0.6, 0.6, 7)
    ew_grid = 55.0 * 10 ** (0.62 * a_grid)  # maa
    inv = cog_invert(a_grid, ew_grid, ew_meas=50.0, ew_err=2.5)
    out["cog_a"] = inv["a"]
    out["cog_sigma_a"] = inv["sigma_a"]
    out["cog_dadew"] = inv["dadew"]
    out["cog_out_of_range"] = float(inv["out_of_range"])

    agg = aggregate(AGG_LINES)
    for key, val in agg.items():
        out[f"agg_{key}"] = float(val)
    single = aggregate(AGG_LINES[:1])
    out["agg_single_a_mean"] = single["a_mean"]
    out["agg_single_sigma_mean"] = single["sigma_mean"]

    out["rew_50_5000"] = reduced_ew(50.0, 5000.0)
    out["rew_gui_50_5000"] = gui_ab.reduced_ew(50.0, 5000.0)  # must be the same function

    # GUI error propagation + velocity screen (move to tumult.abundance in Phase 2)
    x, y = gauss_line()
    prof = gui_ab.measure_ew(Spectrum(x, y), window=(4999.0, 5001.0), snr=80,
                             cont_pct=2.0, profile="gauss", x0=LINE_X0)
    out["measure_ew_mAA"] = prof.maa
    out["measure_ew_err_phot_mAA"] = prof.err_phot_maa
    out["measure_ew_err_cont_mAA"] = prof.err_cont_maa
    sa, sap, sac = gui_ab.sigma_a_from_rew(prof, drewda=0.4)
    out["sigma_a"] = sa
    out["sigma_a_phot"] = sap
    out["sigma_a_cont"] = sac
    out["velocity_hwhm"] = gui_ab.velocity_hwhm(0.22, 5000.0)
    return out


# --- dispatch run reading (mu-integrated flux, mode detection) --------------------------------

DEMO_RUN = ROOT / "data" / "demo"


def run_reading_cases():
    import tumult

    run = tumult.read(str(DEMO_RUN))
    out = {
        "lam": np.asarray(run.lam, dtype=np.float64),
        "flux": np.asarray(run.flux, dtype=np.float64),
        "flux_lte": np.asarray(run.flux_lte, dtype=np.float64),
        "flux_cnt": np.asarray(run.flux_cnt, dtype=np.float64),
        "wts": np.asarray(run.wts, dtype=np.float64),
        "mus": np.asarray(run.mus, dtype=np.float64),
        "ie": np.asarray(run.ie, dtype=np.float64),
        "mode": np.array(str(run.mode)),
        "legacy": np.array(float(bool(run.legacy))),
    }
    spec = run.get_spectrum()
    out["spectrum_x"] = np.asarray(spec.xdata.value, dtype=np.float64)
    out["spectrum_y"] = np.asarray(getattr(spec.ydata, "value", spec.ydata), dtype=np.float64)
    spec_lte = run.get_spectrum(imu=0, LTE=True)
    out["spectrum_lte_mu0_y"] = np.asarray(
        getattr(spec_lte.ydata, "value", spec_lte.ydata), dtype=np.float64)
    return out


# --- model-atom reading (AtomClass slow reader + FastAtom pandas reader) ----------------------

DEMO_ATOM = ROOT / "data" / "demo" / "atoms" / "atom.O"


def atom_cases():
    from tumult.atom import AtomClass
    from tumult.fast_atom import FastAtom

    atom = AtomClass(str(DEMO_ATOM), read_bf=True)
    out = {
        "nk": float(atom.nk),
        "nline": float(atom.nline),
        "ncont": float(atom.ncont),
        "abnd": atom.abnd,
        "awgt": atom.awgt,
        "level_ev": np.array([lvl.ev for lvl in atom.level]),
        "level_g": np.array([lvl.g for lvl in atom.level]),
        "line_lam0": np.array([t.lam0 for t in atom.trans[: atom.nline]]),
        "line_bij": np.array([t.bij for t in atom.trans[: atom.nline]]),
        # partition function leans on EK and saha_const - the constants-consolidation canary
        "pf_5000_1e12": atom.partitionFunction(5000.0, 1e12),
        "pf_8000_1e13": atom.partitionFunction(8000.0, 1e13),
    }
    fast = FastAtom(str(DEMO_ATOM))
    out["fast_lam0"] = np.asarray(fast.bbt.lam0.values, dtype=float)
    out["fast_ev"] = np.asarray(fast.lvl.ev.values, dtype=float)
    return out


# --- standalone HOP library (EOS + opacity through libhop.so) --------------------------------

HOP_LIB = ROOT / "HOP" / "standalone" / "build" / "libhop.so"
HOP_LINELIST = ROOT / "input" / "linelists" / "master_line.txt"

# a coarse T x rho node grid (4 x 6 = 24 cells) and a narrow 10 AA window: small enough to
# solve in seconds, wide enough to cross the Mg b lines and a good deal of Fe I
HOP_DLAM = 10.0        # line-wing search radius [AA]
HOP_VMIC = 1.0         # [km/s]; must be pinned, it moves the line opacity a lot


def hop_grid():
    """The (temp, rho) node grid the HOP cases solve on, both (ntemp, nrho)."""
    from tumult.hop import logstep

    tab_temp = logstep(3500.0, 12000.0, 0.35)
    tab_rho = logstep(1e-10, 1e-7, 1.2)
    return np.meshgrid(tab_temp, tab_rho, indexing="ij")


def hop_wavelengths():
    """The vacuum wavelength grid the opacity cases use."""
    from tumult.hop import wavelength_grid

    return wavelength_grid(5160.0, 5170.0, R=20000)


def hop_eos_cases():
    from tumult.hop import HopEOS

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


def hop_opacity_cases():
    from tumult.hop import HopEOS

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


ALL_CASES = {
    "ew_profiles": ew_profile_cases,
    "convolutions": convolution_cases,
    "airvac": airvac_cases,
    "abundance": abundance_cases,
    "run_reading": run_reading_cases,
    "atom": atom_cases,
    "hop_eos": hop_eos_cases,
    "hop_opacity": hop_opacity_cases,
}
