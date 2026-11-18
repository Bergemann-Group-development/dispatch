"""Figure of the mesh grading (Sec. "Mesh"): how amr='P' places the depth points of 3D
snapshots, and the temperature against log tau_500 on the native and graded mesh, one row
per snapshot.

Each row reads two 3D LTE runs of one M3DIS snapshot, one on the snapshot's own mesh
(nx=ny=nz=-1, amr='') and one graded (nx=ny=-1, nz=96, amr='P'), both with dims=32. The
solar runs need about 40 GiB each, so run them on a cluster node and this script there,
from the experiment directory:

    python docs/manual/figs/mesh_grading.py NATIVE_RUN AMR_RUN [NATIVE_RUN AMR_RUN ...]
"""
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import tumult as tm
from science_case import FIGS, wide


def accumulated_lpe(atmos):
    """Depth below log tau_500 = 0 and the accumulated |change| of the horizontally averaged
    ln Pe from the top, normalised to 1."""
    lpe = np.log(atmos.pe).mean(axis=(0, 1))[::-1]
    acc = np.concatenate([[0], np.cumsum(np.abs(np.diff(lpe)))])
    return -atmos.height[::-1], acc / acc[-1]


def stellar_params(run):
    """T_eff, log g and [Fe/H] as the snapshot's file name gives them."""
    name = Path(run.nml["atmos_params"]["atmos_file"]).name
    teff, logg, feh = re.search(r"_t(\d+)_g([+-]?[\d.]+)_z([+-]?[\d.]+)", name).groups()
    return "\n".join((rf"$T_{{\rm eff}} = {int(teff)}\,$K", rf"$\log_{{10}} g = {float(logg):.2f}$",
                      rf"[Fe/H]$\,= {float(feh):+.2f}$"))


def main(pairs):
    fig, rows = plt.subplots(len(pairs), 3, figsize=wide(len(pairs)), layout="constrained",
                             squeeze=False)
    for axes, (native, graded) in zip(rows, pairs):
        plot_row(axes, native, graded)
    for letter, ax in zip("abcdefghijkl", rows.flat):
        ax.annotate(f"({letter})", xy=(0, 1), xycoords="axes fraction", xytext=(8, -8),
                    textcoords="offset points", ha="left", va="top")
    fig.savefig(FIGS / "mesh_grading.pdf")


def plot_row(axes, native, graded):
    ax1, ax2, ax3 = axes
    ax3.sharex(ax2)
    ax3.sharey(ax2)

    depth, acc = accumulated_lpe(native.atmos)
    dims = native.nml["atmos_params"]["dims"]
    for level in np.linspace(0, 1, dims + 1):
        ax1.axhline(level, color="lightgrey", lw=0.4)
    ax1.plot(depth, acc, color="black")
    points = -graded.atmos.height
    ax1.plot(points, np.interp(points, depth, acc), ".", color="black")
    ax1.set_xlabel("Depth [Mm]")
    ax1.set_ylabel(r"$\sum|\Delta\langle\ln P_{\rm e}\rangle|$, normalised")
    # the curve rises from the lower left to the upper right, leaving the lower right free
    ax1.annotate(stellar_params(native), xy=(1, 0), xycoords="axes fraction", xytext=(-8, 8),
                 textcoords="offset points", ha="right", va="bottom", ma="right")

    # both panels share the range of the native mesh, which is the larger one
    ltau, temp = native.atmos.ltau, native.atmos.temp
    lims = dict(xmin=ltau.min(), xmax=ltau.max(), ymin=temp.min(), ymax=temp.max())
    edges = np.linspace(ltau.min(), ltau.max(), 61)
    hists = []
    for ax, run, label in ((ax2, native, "native"), (ax3, graded, "amr='P'")):
        atmos = run.atmos
        plt.sca(ax)
        atmos.plot_ltau_cube("temp", density=True, cmap="GnBu", bins=300, **lims)
        atmos.plot_ltau_mean("temp", color="black")
        ax.annotate(f"{label}, {atmos.ltau.shape[-1]} points", xy=(0, 1),
                    xycoords="axes fraction", xytext=(8, -19), textcoords="offset points",
                    ha="left", va="top")
        ax.set_xlabel(r"$\log_{10}\tau_{500}$")

        # the histogram of the cells in log tau_500 stands on the panel's top axis
        hist = ax.inset_axes((0, 1, 1, 0.25), sharex=ax)
        hists.append(hist)
        hist.hist(atmos.ltau.ravel(), bins=edges, rwidth=0.8, color="black")
        hist.axis("off")
        hist.set_ylim(bottom=0)     # the columns stand on the panel's top axis
    hists[1].sharey(hists[0])
    ax2.set_ylabel(r"$T$ [K]")
    ax3.tick_params(labelleft=False)


if __name__ == "__main__":
    runs = [tm.read(arg) for arg in sys.argv[1:]]
    main(list(zip(runs[::2], runs[1::2])))
