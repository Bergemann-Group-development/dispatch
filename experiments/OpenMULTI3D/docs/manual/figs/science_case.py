"""Figures of the manual drawn with the plotting methods of tumult.

Most show the worked example, Mg I in the solar flux spectrum; the reference sections
reuse its runs for the Grotrian diagram, the contribution functions and the
centre-to-limb variation, and a 3D LTE run of the solar snapshot for the atmosphere.

Run from the experiment directory with the virtual environment active:

    python docs/manual/figs/science_case.py

It launches the 1D NLTE runs of the example (a few seconds each on a laptop) and the 3D
run (about 10 s on eight threads) into data/manual/ and writes the figures next to this
script. The observation is the KPNO solar flux atlas of Kurucz (2005) in input/spectra/ (air
wavelengths, normalised flux; input_data dataset kpno).
"""
from pathlib import Path

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import tumult as tm
from tumult.linemask import parse_linemask_file

FIGS = Path(__file__).resolve().parent
ATLAS = Path("input/spectra/KPNO_4000_9000.txt")
ATMOSPHERE = "./input/atmos/sun.mod"
SNAPSHOT = "input/atmos/SunLTE_t5777_g4.440_z0.000_sn779_m3dis"
ATOM = "./input/atoms/atom.mg86g"
WINDOW_FILE = Path("data/manual/mg_windows.txt")
# synthesis windows, keyed by the line they are centred on
WINDOWS = {"5172.7": (5170.2, 5175.2), "5528.4": (5527.4, 5529.4),
           "5711.1": (5710.1, 5712.1)}
# the abundance lines and the windows their equivalent widths are measured over
MASK = parse_linemask_file("input/linemasks/Mg/mg-lmask.txt")
ABUNDANCES = np.round(np.arange(7.28, 7.89, 0.1), 2)
REFERENCE = 7.58
VROT = 1.8                      # solar rotation [km/s], fixed in the broadening fit
FIT_WINDOW = (5710.6, 5711.6)   # 5711.1, used to fit macroturbulence and shift
BLEND_WINDOW = (5171.2, 5174.2)
# levels whose departure coefficients are shown: those of the three lines
LEVELS = {0: r"3s$^2$ $^1$S", 2: r"3p $^3$P$^\circ_1$", 5: r"4s $^3$S",
          4: r"3p $^1$P$^\circ$", 12: r"5s $^1$S", 13: r"4d $^1$D", 85: r"Mg II 3s"}

# the manual's figures use the house style and its colour cycle
plt.style.use(Path(tm.__file__).parent / "Bergemann2020.mplstyle")
CYCLE = plt.rcParams["axes.prop_cycle"].by_key()["color"]
NLTE, LTE = CYCLE[0], CYCLE[1]
OBSERVED = dict(color="lightgrey", marker=".", ls="none", ms=3)
GUIDE = dict(color="lightgrey", ls="--", lw=0.8, zorder=0)


def launch(name, abundance=REFERENCE, cntrbf=False):
    """One 1D NLTE run of the Mg I windows on the solar MARCS model."""
    WINDOW_FILE.parent.mkdir(parents=True, exist_ok=True)
    WINDOW_FILE.write_text("# blue_AA  red_AA\n"
                           + "".join(f"{b} {r}\n" for b, r in WINDOWS.values()))
    cfg = {
        "io_params": {"datadir": "."},
        "atmos_params": {"atmos_file": ATMOSPHERE},
        "atom_params": {"atom_file": ATOM, "abundance": abundance},
        "opac_table_params": {"absmet_file": "./input/absmet_files/absmet"},
        "linelist_params": {"linelist_file": "./input/linelists/master_line.txt"},
        "spectrum_params": {"window_file": f"./{WINDOW_FILE}", "R_sampling": 5e5},
        "output_params": {"save_cntrbf": cntrbf},
    }
    folder, _ = tm.launch(f"data/manual/{name}", nml_settings=cfg,
                          default=[tm.preset("1d"), tm.preset("nlte")],
                          threads=8, verbose=0)
    return tm.read(folder)


def launch_demo(name, blends=False):
    """The Ba II 4554 run of demo.nml; with blends, the version of Sec. 4 that adds a
    line list, a 4550-4558 A window and the solar metal-opacity table."""
    cfg = {
        "io_params": {"datadir": "."},
        "atmos_params": {"atmos_file": ATMOSPHERE},
        "atom_params": {"atom_file": "./input/atoms/atom.ba06", "use_atom_abnd": True},
        "nlte_params": {"convlim": 1e-2},
    }
    if blends:
        cfg |= {
            "opac_table_params": {"absmet_file": "./input/absmet_files/absmet"},
            "linelist_params": {"linelist_file": "./input/linelists/master_line.txt"},
            "spectrum_params": {"R_sampling": 5e5, "aa_blue": 4550, "aa_red": 4558},
        }
    folder, _ = tm.launch(f"data/manual/{name}", nml_settings=cfg, threads=8, verbose=0)
    return tm.read(folder)


def fig_demo(run, blended):
    """The first run: Ba II 4554 in NLTE and LTE on the atom's grid, and the same line
    in a window with the blends of the line list."""
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(6.6, 2.5),
                                        width_ratios=(1, 1, 1.6), layout="constrained")
    plt.sca(ax0)
    for lvl, label, ls in ((1, r"1: 6s $^2$S", "-"), (4, r"4: 6p $^2$P$^\circ_{3/2}$", "--")):
        run.plot_depart(lvl, label=label, ls=ls)
    ax0.axhline(1, **GUIDE)
    ax0.set_xlim(-5, 1)
    ax0.set_ylabel(r"$b$")
    ax0.legend(loc="lower right", fontsize=7)
    ax0.set_title("departures")

    plt.sca(ax1)
    line = run.get_line(4554.033)
    for lte, colour, ls, label in ((False, NLTE, "-", "NLTE"), (True, LTE, "--", "LTE")):
        line.crop(LTE=lte).plot(color=colour, ls=ls, label=label, set_labels=False)
    ax1.set_xlabel(r"$\lambda$ [$\AA$]")
    ax1.set_ylabel(r"$F/F_\mathrm{c}$")
    ax1.set_xlim(4553.5, 4554.6)
    ax1.set_xticks([4553.7, 4554.3])
    ax1.ticklabel_format(axis="x", useOffset=False)
    ax1.legend(loc="center right", handlelength=1.5, fontsize=7)
    ax1.set_title(r"Ba II 4554 $\AA$")

    plt.sca(ax2)
    blended.crop(window=(4550, 4558)).plot(color=NLTE, set_labels=False)
    blended.plot_line_labels(ntop=6)
    ax2.set_xlabel(r"$\lambda$ [$\AA$]")
    ax2.set_xlim(4550, 4558)
    ax2.ticklabel_format(axis="x", useOffset=False)
    ax2.set_title("with line list and window")
    fig.savefig(FIGS / "demo_ba.pdf")
    return {lte: line.crop(LTE=lte).measure(profile="trapz").maa for lte in (False, True)}


def launch_snapshot():
    """A 3D LTE run of the solar snapshot on two continuum wavelengths, for its atmosphere."""
    cfg = {
        "io_params": {"datadir": "."},
        "atmos_params": {"atmos_format": "m3dis", "atmos_file": SNAPSHOT},
        "spectrum_params": {"aa_blue": 5000.0, "aa_red": 5001.0, "nnu": 2},
    }
    folder, _ = tm.launch("data/manual/sun3d", nml_settings=cfg,
                          default=[tm.preset("3d")], threads=8, verbose=0)
    return tm.read(folder)


def fig_grotrian():
    """Grotrian diagram of Mg I, read from the atom file."""
    atom = tm.read_atom(ATOM)
    fig = plt.figure(figsize=(6.6, 2.9), layout="constrained")
    atom.plot_grotrian(ion=1, color="#c9c8c2", lw=0.4)
    plt.gca().tick_params(axis="x", labelsize=7)
    fig.savefig(FIGS / "mg_grotrian.pdf")


def fig_depart(run):
    """Departure coefficients of the lines' levels."""
    fig = plt.figure(figsize=(6.6, 2.8), layout="constrained")
    # colours of the style cycle repeat, so each level also gets its own dash pattern
    styles = ("-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 3)))
    colours = (CYCLE[0], CYCLE[1], CYCLE[2], CYCLE[3], CYCLE[4], CYCLE[5], "k")
    for (lvl, name), ls, colour in zip(LEVELS.items(), styles, colours):
        run.plot_depart(lvl, label=f"{lvl}: {name}", ls=ls, color=colour)
    ax = plt.gca()
    ax.axhline(1, **GUIDE)
    ax.set_title("")
    ax.set_ylabel(r"$b = n_\mathrm{NLTE}/n_\mathrm{LTE}$")
    ax.set_xlim(-5.5, 1.5)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(FIGS / "mg_depart.pdf")


def fig_lines(run, atlas):
    """Broadening fitted to 5711.1, and the blends of Mg b2 labelled. Returns the
    fitted broadening."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.9), width_ratios=(1, 1.4),
                                   layout="constrained")
    plt.sca(ax1)
    obs = atlas.crop(window=FIT_WINDOW)
    obs.plot(label="KPNO atlas", set_labels=False, **OBSERVED)
    best, fitted = run.fit(obs, fit=["zeta_rt", "vshift"], vrot=VROT)
    fitted.plot(color=NLTE, label="NLTE fit", set_labels=False)
    ax1.set_xlim(*FIT_WINDOW)
    ax1.set_ylim(0.1, 1.08)
    ax1.set_xlabel(r"$\lambda$ [$\AA$]")
    ax1.set_ylabel(r"$F/F_\mathrm{c}$")
    ax1.set_xticks([5710.8, 5711.1, 5711.4])
    ax1.ticklabel_format(axis="x", useOffset=False)
    ax1.legend(loc="lower left")
    ax1.set_title("Mg I 5711.1 $\\AA$, fit", pad=4)

    plt.sca(ax2)
    atlas.crop(window=BLEND_WINDOW).plot(label="KPNO atlas", set_labels=False, **OBSERVED)
    run.crop(window=BLEND_WINDOW, vrot=VROT, **best).plot(color=NLTE, label="NLTE",
                                                          set_labels=False)
    ax2.set_xlim(*BLEND_WINDOW)
    ax2.set_xticks([5172, 5173, 5174])
    ax2.set_xlabel(r"$\lambda$ [$\AA$]")
    ax2.ticklabel_format(axis="x", useOffset=False)
    run.plot_line_labels(show_atom=False, ntop=10)
    fig.savefig(FIGS / "mg_lines.pdf")
    return dict(vrot=VROT, **best)


def fig_spectra(run, atlas, broadening):
    """NLTE and LTE profiles of the three windows with the observed atlas."""
    fig, (axes, diffs) = plt.subplots(2, 3, figsize=(6.6, 3.4), sharex="col",
                                      sharey="row", height_ratios=(3, 1),
                                      layout="constrained")
    for ax, diff, (line, win) in zip(axes, diffs, WINDOWS.items()):
        atlas.crop(window=win).plot(ax=ax, label="KPNO atlas", set_labels=False,
                                    **OBSERVED)
        run.crop(window=win, LTE=True, **broadening).plot(
            ax=ax, color=LTE, ls="--", label="LTE", set_labels=False)
        run.crop(window=win, **broadening).plot(
            ax=ax, color=NLTE, label="NLTE", set_labels=False)
        nlte = run.crop(window=win, **broadening)
        lte = run.crop(window=win, LTE=True, **broadening)
        diff.plot(nlte.xdata.value, nlte.ydata.value - lte.ydata.value, color=NLTE)
        diff.axhline(0, **GUIDE)
        ax.set_title(rf"Mg I {line} $\AA$")
        diff.set_xlabel(r"$\lambda$ [$\AA$]")
        ax.set_xlim(*win)
        diff.ticklabel_format(axis="x", useOffset=False)
        diff.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    axes[0].set_ylabel(r"$F/F_\mathrm{c}$")
    diffs[0].set_ylabel("NLTE$-$LTE")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=3)
    fig.savefig(FIGS / "mg_spectra.pdf")


def fig_cntrbf(run):
    """Contribution function of Mg b2 over wavelength and depth, and at selected
    wavelengths, compared with a weak line."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.9), width_ratios=(1.2, 1),
                                   layout="constrained")
    fig.get_layout_engine().set(w_pad=0.15)
    top = run.atmos.ltau.min()
    plt.sca(ax1)
    mesh = run.plot_cntrbf_map(xmin=5171.2, xmax=5174.2, cmap="Blues", title="",
                               norm=matplotlib.colors.PowerNorm(0.5))
    fig.colorbar(mesh, ax=ax1,
                 label=r"$\mathrm{d}(1-F/F_\mathrm{c})\,/\,\mathrm{d}\log\tau_{500}$")
    ax1.set_ylim(1, top)
    ax1.set_title(r"Mg I 5172.7 $\AA$, depression")

    plt.sca(ax2)
    for lam0, label, colour, ls in ((5172.68, "5172.7 core", NLTE, "-"),
                                    (5172.78, r"5172.7 $+0.1\,\AA$", NLTE, "--"),
                                    (5172.98, r"5172.7 $+0.3\,\AA$", NLTE, "-."),
                                    (5173.28, r"5172.7 $+0.6\,\AA$", NLTE, ":"),
                                    (5711.09, "5711.1 core", LTE, "-")):
        run.plot_cntrbf(lam0=lam0, norm=True, color=colour, ls=ls, label=label, title="")
    ax2.set_ylabel("Normalised contribution")
    ax2.set_xlim(top, 1)
    ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax2.set_title("Formation depths")
    fig.savefig(FIGS / "mg_cntrbf.pdf")


def fig_cog(runs, atlas):
    """NLTE and LTE curves of growth with the observed equivalent widths."""
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), layout="constrained")
    results = {}
    for ax, line in zip(axes, MASK):
        win = line.fit_window
        nlte = [r.crop(window=win).measure(profile="trapz").maa for r in runs]
        lte = [r.crop(window=win, LTE=True).measure(profile="trapz").maa for r in runs]
        ew_obs = atlas.measure(win, profile="trapz").maa
        abnd_nlte = tm.cog_invert(ABUNDANCES, nlte, ew_obs)["a"]
        abnd_lte = tm.cog_invert(ABUNDANCES, lte, ew_obs)["a"]
        results[line.lam0] = (ew_obs, abnd_lte, abnd_nlte)

        ax.plot(ABUNDANCES, nlte, "o-", color=NLTE, ms=4, label="NLTE")
        ax.plot(ABUNDANCES, lte, "s--", color=LTE, ms=4, label="LTE")
        ax.axhline(ew_obs, **GUIDE, label="KPNO atlas")
        for a in (abnd_nlte, abnd_lte):
            ax.axvline(a, **GUIDE)
        ax.set_title(rf"Mg I {line.lam0:.1f} $\AA$")
        ax.set_xlabel(r"$A(\mathrm{Mg})$")
        ax.set_ylabel(r"$W$ [m$\AA$]")
    axes[0].legend(loc="upper left")
    fig.savefig(FIGS / "mg_cog.pdf")
    return results


def fig_cog_corr(run, runs):
    """NLTE correction of 5528.4 on the LTE curve of growth of the abundance ladder.
    Returns the correction."""
    line = next(iter(MASK))
    fig = plt.figure(figsize=(3.3, 2.6), layout="constrained")
    run.plot_cog_corr("Mg", LTE_runs=runs, window=line.fit_window)
    corr = run.get_cog_corr("Mg", runs, window=line.fit_window)[0]
    plt.title(rf"Mg I {line.lam0:.1f} $\AA$")
    plt.legend(loc="upper left")
    fig.savefig(FIGS / "mg_cog_corr.pdf")
    return corr


def fig_clv(run):
    """Centre-to-limb variation of the equivalent width of 5711.1."""
    fig = plt.figure(figsize=(3.3, 2.6), layout="constrained")
    line = run.get_line(5711.09)
    line.plot_clv(dlam=0.3, LTE=True, title="", color=LTE, marker="s", mfc="none", ms=7,
                  label="LTE")
    line.plot_clv(dlam=0.3, title="", color=NLTE, marker="o", ms=4, label="NLTE")
    plt.legend(loc="lower right")
    fig.savefig(FIGS / "mg_clv.pdf")
    return line.get_clv(dlam=0.3)


def fig_snapshot(run):
    """Temperature of the solar 3D snapshot: horizontal and vertical slices, against
    log tau_500, and its distribution at log tau_500 = 0."""
    atmos = run.atmos
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 5.0), layout="constrained")
    plt.sca(axes[0, 0])
    mesh = atmos.plot_slice("temp", iz=int(np.argmin(np.abs(atmos.height))),
                            cmap="inferno", title=r"$z=0$")
    fig.colorbar(mesh, ax=axes[0, 0], label=r"$T$ [K]")
    axes[0, 0].set_aspect("equal")

    plt.sca(axes[0, 1])
    mesh = atmos.plot_slice("temp", iy=0, cmap="inferno", title=r"$y=0$")
    fig.colorbar(mesh, ax=axes[0, 1], label=r"$T$ [K]")

    plt.sca(axes[1, 0])
    atmos.plot_ltau_cube("temp", color="#c9c8c2")
    atmos.plot_ltau_mean("temp", color=NLTE, label=r"mean on $\tau_{500}$ surfaces")
    axes[1, 0].set_xlim(-5, 3)
    axes[1, 0].set_ylim(3000, 14000)
    axes[1, 0].set_xlabel(r"$\log\tau_{500}$")
    axes[1, 0].set_ylabel(r"$T$ [K]")
    axes[1, 0].legend(loc="lower right")

    plt.sca(axes[1, 1])
    atmos.plot_ltau_hist(0, "temp", show_mean=True, color=NLTE, bins=30)
    axes[1, 1].set_xlabel(r"$T$ at $\log\tau_{500}=0$ [K]")
    axes[1, 1].set_ylabel("Columns")
    fig.savefig(FIGS / "sun3d_atmos.pdf")


def main():
    demo = launch_demo("demo")
    print("demo W [mA]:", fig_demo(demo, launch_demo("demo_blends", blends=True)))

    obs = tm.load_spectrum(ATLAS)
    reference = launch("mg_sun", cntrbf=True)
    ladder = [launch(f"mg_sun_A{round(a*100)}", abundance=a) for a in ABUNDANCES]

    fig_grotrian()
    fig_depart(reference)
    broadening = fig_lines(reference, obs)
    fig_spectra(reference, obs, broadening)
    fig_cntrbf(reference)
    results = fig_cog(ladder, obs)
    corr = fig_cog_corr(reference, ladder)
    mus, clv = fig_clv(reference)
    fig_snapshot(launch_snapshot())

    print("broadening:", broadening)
    win = WINDOWS["5172.7"]
    for kind, lte in (("NLTE", False), ("LTE", True)):
        core = reference.crop(window=win, LTE=lte, **broadening).ydata.value.min()
        print(f"5172.7 core depth {kind}: {core:.3f}")
    print(f"5172.7 core observed: {obs.crop(window=(5172.5, 5172.9)).ydata.value.min():.3f}")
    for line in MASK:
        nlte = reference.crop(window=line.fit_window).measure(profile="trapz").maa
        lte = reference.crop(window=line.fit_window, LTE=True).measure(profile="trapz").maa
        print(f"{line.lam0:.1f}: W_NLTE={nlte:6.1f} mA  W_LTE={lte:6.1f} mA")
    for lam0, (ew_obs, abnd_lte, abnd_nlte) in results.items():
        print(f"{lam0:.1f}: W_obs={ew_obs:6.1f} mA  A_LTE={abnd_lte:.2f}  A_NLTE={abnd_nlte:.2f}"
              f"  correction={abnd_nlte - abnd_lte:+.2f}")
    print(f"5528.4 plot_cog_corr correction at A=7.58: {corr:+.3f}")
    print("CLV 5711.1:", dict(zip(np.round(mus, 3), np.round(clv, 1))))


if __name__ == "__main__":
    main()
