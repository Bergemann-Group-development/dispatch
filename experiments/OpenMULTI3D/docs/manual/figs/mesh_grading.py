"""Figure of the mesh grading (Sec. "Mesh"): how amr='P' places the depth points of the
solar 3D snapshot, and the temperature against log tau_500 on the native and graded mesh.

It reads two 3D LTE runs of the snapshot of science_case.py with save_atmos=T, one on the
snapshot's own mesh (nx=ny=nz=-1, amr='') and one graded (nx=ny=-1, nz=96, amr='P'),
both with dims=32. They need about 40 GiB each, so run them on a cluster node and this
script there, from the experiment directory:

    python docs/manual/figs/mesh_grading.py NATIVE_RUN AMR_RUN
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

import tumult as tm
from science_case import FIGS, wide

GREY = "#c9c8c2"


def accumulated_lpe(atmos):
    """Accumulated |change| of the horizontally averaged ln Pe from the top, normalised to 1."""
    lpe = np.log(atmos.pe).mean(axis=(0, 1))[::-1]
    acc = np.concatenate([[0], np.cumsum(np.abs(np.diff(lpe)))])
    return atmos.height[::-1], acc / acc[-1]


def main(native, graded):
    fig, axes = plt.subplots(1, 3, figsize=wide(), layout="constrained")
    ax1, ax2, ax3 = axes
    ax3.sharex(ax2)
    ax3.sharey(ax2)

    height, acc = accumulated_lpe(native.atmos)
    dims = native.nml["atmos_params"]["dims"]
    for level in np.linspace(0, 1, dims + 1):
        ax1.axhline(level, color="lightgrey", lw=0.4)
    ax1.plot(height, acc, color="black")
    points = graded.atmos.height
    ax1.plot(points, np.interp(-points, -height, acc), ".", color="black")
    ax1.set_xlabel("Height [Mm]")
    ax1.set_ylabel(r"$\sum|\Delta\langle\ln P_{\rm e}\rangle|$, normalised")
    ax1.invert_yaxis()

    for ax, run, label in ((ax2, native, "native mesh"), (ax3, graded, "amr='P'")):
        atmos = run.atmos
        plt.sca(ax)
        atmos.plot_ltau_cube("temp", color=GREY)
        ax.plot(atmos.ltau[0, 0], atmos.temp[0, 0], ".", color="black")
        ax.set_title(f"{label}, {atmos.ltau.shape[-1]} depth points")
        ax.set_xlabel(r"$\log_{10}\tau_{500}$")
    ax2.autoscale()     # plot_ltau_cube fixes the limits to each run's own range
    ax2.set_ylabel(r"$T$ [K]")
    ax3.tick_params(labelleft=False)

    for letter, ax in zip("abc", axes):
        ax.annotate(f"({letter})", xy=(0, 1), xycoords="axes fraction", xytext=(8, -8),
                    textcoords="offset points", ha="left", va="top")
    fig.savefig(FIGS / "mesh_grading.pdf")


if __name__ == "__main__":
    main(tm.read(sys.argv[1]), tm.read(sys.argv[2]))
