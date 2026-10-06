"""The campaign's <3D> average of one Stagger model, from its full-resolution snapshots.

    python avg3d.py OUT.atmos FEH ALPHA SNAPSHOT [SNAPSHOT ...] [--real Q] [--top TOP]

The snapshots (Stagger-grid-2 atm3d.* files, mesh.* beside them) are read with
tumult.hop.read_atmosphere, their columns pooled and averaged with tumult.avg3d
(reports/2026-10-06_avg3d_study.md): surfaces of log tau_500, T = <T^4>^(1/4), geometric
rho and P_e, and the surface's log tau_500 as the depth column.  tau_500 comes from HOP's
continuum opacity at [Fe/H] = FEH and [alpha/Fe] = ALPHA, tabulated in (T, rho) once per
model.  No dispatch.x run and no resampling are involved.

The top surface is the one that a fraction Q of the columns reach (default 0.5, the median
column top, at most log tau_500 = -5).  Shorter columns are continued above their top as
TOP says (tumult.avg3d._on_surfaces): 'isothermal' (default) at their top temperature,
'gradient' along the mean T, rho and P_e of the complete columns
(reports/2026-10-06_phase7.md compares them).

The result is dispatch.x's 'text' format (log tau_500, T, P_e, rho, vmic), read with
atmos_format='text' use_ne=F amr='T': rho given, n_e from the EOS.  vmic is the rms of
the three velocity components on the surface divided by sqrt(3), the convention of the
Stagger MULTI-format <3D> files.  Each call appends one summary line to OUT.atmos.txt:
the settings, the top and bottom surfaces, the share of extended columns at the top
and the snapshots.
"""
import argparse
from pathlib import Path

from tumult import avg3d
from tumult.hop import read_atmosphere

RECIPE = dict(surface='tau500', recipe='t4')


def average(snapshots, feh, alpha, real=0.5, top='isothermal'):
    atm = [read_atmosphere(str(f), 'stagger2') for f in snapshots]
    table = avg3d.OpacityTable.for_snapshots(atm, feh=feh, alpha_enhancement=alpha)
    cols = avg3d.Columns.concatenate([avg3d.columns(a, table) for a in atm])
    return avg3d.surface_average(cols, real=real, top=top, **RECIPE)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('out', type=Path)
    p.add_argument('feh', type=float)
    p.add_argument('alpha', type=float)
    p.add_argument('snapshots', nargs='+', type=Path)
    p.add_argument('--real', type=float, default=0.5)
    p.add_argument('--top', choices=avg3d.TOPS, default='isothermal')
    a = p.parse_args(argv)
    av = average(a.snapshots, a.feh, a.alpha, a.real, a.top)
    avg3d.write_text(a.out, av, a.feh, depth='tau')
    line = (f'{a.out.name} snapshots={len(a.snapshots)} real={a.real} top={a.top} '
            f'ltau_top={av["x"][0]:.3f} ltau_bottom={av["x"][-1]:.3f} '
            f'extended_at_top={av["extrapolated"][0]:.3f} '
            f'{" ".join(f.name for f in a.snapshots)}')
    with open(a.out.with_suffix(a.out.suffix + '.txt'), 'a') as f:
        f.write(line + '\n')
    print(line)


if __name__ == '__main__':
    main()
