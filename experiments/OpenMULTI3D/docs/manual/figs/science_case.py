"""Figures of the manual's worked example: Mg I in the solar flux spectrum.

Run from the experiment directory with the virtual environment active:

    python docs/manual/figs/science_case.py [--atlas KPNO_4000_9000.txt]

It launches the 1D NLTE runs of the example (a few seconds each on a laptop) into
data/manual/ and writes the figures next to this script. ``--atlas`` names a
two-column solar flux atlas (wavelength in air, normalised flux), such as the KPNO
atlas of Kurucz (2005); without it the spectra and curves of growth are drawn without
the observation.
"""
import argparse
from pathlib import Path

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import tumult as tm

FIGS = Path(__file__).resolve().parent
ATMOSPHERE = ("./input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00"
              "_o+0.00_r+0.00_s+0.00.mod")
ATOM = "./input/atoms/atom.mg86g"
WINDOW_FILE = Path("data/manual/mg_windows.txt")
# synthesis windows, keyed by the line they are centred on
WINDOWS = {5172: (5170.2, 5175.2), 5528: (5527.4, 5529.4), 5711: (5710.1, 5712.1)}
# the abundance lines of input/linemasks/Mg/mg-lmask.txt: centre and the
# window over which the equivalent width is integrated
EW_WINDOWS = {5528: (5528.405, (5528.1, 5528.7)), 5711: (5711.088, (5710.8, 5711.4))}
ABUNDANCES = np.round(np.arange(7.28, 7.89, 0.1), 2)
REFERENCE = 7.58
BROADENING = dict(vrot=1.8, zeta_rt=3.5)   # solar rotation and macroturbulence [km/s]
# levels whose departure coefficients are shown: those of the three lines
LEVELS = {0: r"3s$^2$ $^1$S", 2: r"3p $^3$P$^\circ_1$", 5: r"4s $^3$S",
          4: r"3p $^1$P$^\circ$", 12: r"5s $^1$S", 13: r"4d $^1$D", 85: r"Mg II 3s"}

# the manual's figures use the house style and its colour cycle
plt.style.use(Path(tm.__file__).parent / "Bergemann2020.mplstyle")
CYCLE = plt.rcParams["axes.prop_cycle"].by_key()["color"]
NLTE, LTE = CYCLE[0], CYCLE[1]
INK = "black"


def launch(name, abundance=REFERENCE, cntrbf=False):
    """One 1D NLTE run of the Mg I windows on the solar MARCS model."""
    WINDOW_FILE.parent.mkdir(parents=True, exist_ok=True)
    WINDOW_FILE.write_text("# blue_AA  red_AA\n"
                           + "".join(f"{b} {r}\n" for b, r in WINDOWS.values()))
    settings = {
        "io_params": {"datadir": "."},
        "atmos_params": {"atmos_file": ATMOSPHERE},
        "atom_params": {"atom_file": ATOM, "abundance": abundance},
        "opac_table_params": {"absmet_file": "./input/absmet_files/absmet"},
        "linelist_params": {"linelist_file": "./input/linelists/master_line.txt"},
        "spectrum_params": {"window_file": f"./{WINDOW_FILE}", "R_sampling": 5e5},
        "output_params": {"save_cntrbf": cntrbf},
    }
    folder, _ = tm.launch(f"data/manual/{name}", nml_settings=settings,
                          default=[tm.preset("1d"), tm.preset("nlte")],
                          threads=8, verbose=0)
    return tm.read(folder)


def ew(spectrum, line):
    """Equivalent width in mA of a normalised spectrum over the line's mask window."""
    return spectrum.crop(window=EW_WINDOWS[line][1]).get_weq(target_unit=u.mAA).value


def fig_atom(run):
    """Grotrian diagram of Mg I and the departure coefficients of the lines' levels."""
    atom = run.atoms[0]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.6, 5.6), height_ratios=(1, 1),
                                   layout="constrained")
    plt.sca(ax1)
    atom.plot_grotrian(ion=1, color="#c9c8c2", lw=0.4)
    ax1.tick_params(axis="x", labelsize=7)
    ax1.set_title("Mg I levels and transitions")
    plt.sca(ax2)
    for level, name in LEVELS.items():
        run.plot_depart(level, label=name)
    ax2.axhline(1, color="lightgrey", lw=0.6, zorder=0)
    ax2.set_title("Departure coefficients")
    ax2.set_ylabel(r"$b = n_\mathrm{NLTE}/n_\mathrm{LTE}$")
    ax2.set_xlim(1.5, -5.5)
    ax2.legend(loc="lower left", ncol=3)
    fig.savefig(FIGS / "mg_atom.pdf")


def fig_spectra(run, atlas):
    """NLTE and LTE profiles of the three windows with the observed atlas."""
    fig, axes = plt.subplots(1, 3, figsize=(6.6, 2.6), sharey=True, layout="constrained")
    for ax, (line, window) in zip(axes, WINDOWS.items()):
        if atlas is not None:
            atlas.crop(window=window).plot(ax=ax, color=INK, lw=0.8, label="KPNO atlas",
                                           set_labels=False)
        run.crop(window=window, LTE=True, **BROADENING).plot(
            ax=ax, color=LTE, ls="--", label="LTE", set_labels=False)
        run.crop(window=window, **BROADENING).plot(
            ax=ax, color=NLTE, label="NLTE", set_labels=False)
        ax.set_title(rf"Mg I {line} $\AA$")
        ax.set_xlabel(r"$\lambda$ [$\AA$]")
        ax.set_xlim(*window)
        ax.ticklabel_format(axis="x", useOffset=False)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    axes[0].set_ylabel(r"$F/F_\mathrm{c}$")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=3)
    fig.savefig(FIGS / "mg_spectra.pdf")


def fig_cntrbf(run):
    """Contribution function of Mg b2 over wavelength and depth, and at selected
    wavelengths, compared with a weak line."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.9), width_ratios=(1.2, 1),
                                   layout="constrained")
    top = run.atmos.ltau.min()
    plt.sca(ax1)
    mesh = run.plot_cntrbf_map(xmin=5171.2, xmax=5174.2, cmap="Blues", title="",
                               norm=matplotlib.colors.PowerNorm(0.5))
    fig.colorbar(mesh, ax=ax1,
                 label=r"$\mathrm{d}(1-F/F_\mathrm{c})\,/\,\mathrm{d}\log\tau_{500}$")
    ax1.set_ylim(1, top)
    ax1.set_title("Mg I 5172, depression")

    plt.sca(ax2)
    for lam0, label, colour, ls in ((5172.68, "5172 core", NLTE, "-"),
                                    (5172.78, r"5172 $+0.1\,\AA$", NLTE, "--"),
                                    (5172.98, r"5172 $+0.3\,\AA$", NLTE, "-."),
                                    (5173.28, r"5172 $+0.6\,\AA$", NLTE, ":"),
                                    (5711.09, "5711 core", LTE, "-")):
        run.plot_cntrbf(lam0=lam0, norm=True, color=colour, ls=ls, label=label, title="")
    ax2.set_ylabel("Normalised contribution")
    ax2.set_xlim(1, top)
    ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax2.set_title("Formation depths")
    fig.savefig(FIGS / "mg_cntrbf.pdf")


def fig_cog(runs, atlas):
    """NLTE and LTE curves of growth with the observed equivalent widths."""
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), layout="constrained")
    results = {}
    for ax, line in zip(axes, EW_WINDOWS):
        nlte = np.array([ew(run.crop(window=WINDOWS[line]), line) for run in runs])
        lte = np.array([ew(run.crop(window=WINDOWS[line], LTE=True), line) for run in runs])
        ax.plot(ABUNDANCES, nlte, "o-", color=NLTE, ms=4, label="NLTE")
        ax.plot(ABUNDANCES, lte, "s--", color=LTE, ms=4, label="LTE")
        if atlas is not None:
            observed = ew(atlas, line)
            a_nlte = np.interp(observed, nlte, ABUNDANCES)
            a_lte = np.interp(observed, lte, ABUNDANCES)
            results[line] = (observed, a_lte, a_nlte)
            ax.axhline(observed, color=INK, lw=0.8, label="KPNO atlas")
            bottom = ax.get_ylim()[0]
            for a, colour in ((a_nlte, NLTE), (a_lte, LTE)):
                ax.plot([a, a], [bottom, observed], color=colour, lw=0.8)
            ax.set_ylim(bottom=bottom)
        ax.set_title(rf"Mg I {line} $\AA$")
        ax.set_xlabel(r"$A(\mathrm{Mg})$")
        ax.set_ylabel(r"$W$ [m$\AA$]")
    axes[0].legend(loc="upper left")
    fig.savefig(FIGS / "mg_cog.pdf")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--atlas", type=Path, help="solar flux atlas, two columns")
    args = parser.parse_args()
    atlas = tm.load_spectrum(args.atlas) if args.atlas else None

    reference = launch("mg_sun", cntrbf=True)
    ladder = [launch(f"mg_sun_A{round(a*100)}", abundance=a) for a in ABUNDANCES]

    fig_atom(reference)
    fig_spectra(reference, atlas)
    fig_cntrbf(reference)
    results = fig_cog(ladder, atlas)

    for line in EW_WINDOWS:
        print(f"{line}: W_NLTE={ew(reference.crop(window=WINDOWS[line]), line):6.1f} mA"
              f"  W_LTE={ew(reference.crop(window=WINDOWS[line], LTE=True), line):6.1f} mA")
    for line, (observed, a_lte, a_nlte) in results.items():
        print(f"{line}: W_obs={observed:6.1f} mA  A_LTE={a_lte:.2f}  A_NLTE={a_nlte:.2f}"
              f"  correction={a_nlte - a_lte:+.2f}")


if __name__ == "__main__":
    main()
