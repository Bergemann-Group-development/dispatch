#!/usr/bin/env python3
"""Namelists of the <3D> averaging study (reports/2026-10-06_avg3d_study.md).

Every run computes the emergent flux on the MARCS flux sampling (flx_wavelengths.vac,
100724 vacuum wavelengths from 130 nm to 20 micron) with the campaign's FGK line list and
no output scattering iterations, so that the flux integrates to Teff
(tumult.utils.calc_teff).

    python3 avg3d_study.py OUT 3d   # 3D LTE, 3 snapshots of the 4 reference models
    python3 avg3d_study.py OUT 1d   # 1D LTE on every average and reference model

The 1d mode takes the averages that avg3d_study_avg.py wrote to
OUT/avg3d_study/avg/<model>/*.atmos and adds three references per model: the MARCS twin
(marcs), the Stagger <3D> average (stagger, av3d/atmos.<model>) and the campaign's own
average (campaign, avg3d/<model>.atmos).

OUT is a local folder that is copied to $C/work/avg3d_study/; all paths inside the
namelists are relative to $C/work, as in make_nml.py.  Standard library only.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from make_nml import group, linelist_group, namelist, PILOT, HALTON_10  # noqa: E402

LAM_FILE = 'avg3d_study/flx_wavelengths.vac'
SPECTRUM = group('spectrum_params', lam_file=LAM_FILE, in_air=False)
QUIET = dict(save_linelist=False)

# Reference stars: Sun, metal-poor dwarf, giant, metal-poor giant (pilot models).  The
# snapshots are the first three of the campaign's Halton order, so the first two are
# the pilot's 3D snapshots.
MODELS = {
    't5777g44m00': ('t5777/t5777g44m00', ['00459', '00091', '00581']),
    't60g40m20': ('t60/t60g40m20', ['00091', '00041', '00113']),
    't45g25m00': ('t45/t45g25m00_new', ['00206', '00104', '00273']),
    't45g25m20': ('t45/t45g25m20_new', ['00233', '00045', '00287']),
}
PILOT_ROW = {Path(f).name.replace('_new', ''): (marcs, feh, alpha)
             for f, marcs, _, feh, alpha in PILOT}
assert HALTON_10[:3] == [7, 2, 9]


def composition(model):
    _, feh, alpha = PILOT_ROW[model]
    return feh, (dict(alpha_enhancement=alpha) if alpha else None)


def snapshot_file(model, tag):
    folder, _ = MODELS[model]
    stem = folder.split('/')[1].replace('_new', '')
    return f'stagger2/{folder}/atm3d.{stem}05_{tag}_80x80x240'


def flux_3d(model, tag):
    feh, comp = composition(model)
    atmos = dict(atmos_format='stagger2', atmos_file=snapshot_file(model, tag))
    return namelist(atmos, geometry='3d', feh=feh, linelist=None, rt_output=QUIET,
                    composition=comp, extra=[linelist_group('fgk'), SPECTRUM])


def flux_1d(atmos, feh, comp):
    return namelist(atmos, feh=feh, linelist=None, rt_output=QUIET, composition=comp,
                    extra=[linelist_group('fgk'), SPECTRUM])


def references(model):
    return {'marcs': dict(atmos_format='Marcs', atmos_file=PILOT_ROW[model][0]),
            'stagger': dict(atmos_file=f'av3d/atmos.{model}'),
            'campaign': dict(atmos_format='text', atmos_file=f'avg3d/{model}.atmos',
                             use_ne=False, amr='T')}


def write_1d(out):
    tasks = []
    for model in MODELS:
        feh, comp = composition(model)
        runs = references(model)
        for f in sorted((out / 'avg3d_study/avg' / model).glob('*.atmos')):
            runs[f.stem] = dict(atmos_format='text', use_ne=False, amr='T',
                                atmos_file=str(f.relative_to(out)))
        for name, atmos in runs.items():
            rel = Path('avg3d_study/1d') / model / f'{name}.nml'
            (out / rel).parent.mkdir(parents=True, exist_ok=True)
            (out / rel).write_text(flux_1d(atmos, feh, comp))
            tasks.append(str(rel))
    (out / 'avg3d_study/tasks_1d.txt').write_text('\n'.join(tasks) + '\n')
    print(f'{len(tasks)} 1D namelists')


def write_3d(out):
    tasks = []
    for rank in range(3):                      # every model's first snapshot first
        for model, (_, tags) in MODELS.items():
            rel = Path('avg3d_study/3d') / model / f's{tags[rank]}.nml'
            (out / rel).parent.mkdir(parents=True, exist_ok=True)
            (out / rel).write_text(flux_3d(model, tags[rank]))
            tasks.append(str(rel))
    (out / 'avg3d_study/tasks_3d.txt').write_text('\n'.join(tasks) + '\n')
    print(f'{len(tasks)} 3D namelists')


if __name__ == '__main__':
    out = Path(sys.argv[1])
    {'3d': write_3d, '1d': write_1d}[sys.argv[2]](out)
