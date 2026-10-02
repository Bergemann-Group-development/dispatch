#!/usr/bin/env python3
"""Write the dispatch.x namelists of the alpha campaign.

Runs on the cluster login node or locally (standard library only).  Every namelist is
written relative to the campaign work directory, which holds the short-path symlinks the
namelists use:

    input    -> <clone>/input                          (atoms, abundances, line lists)
    marcs    -> storm/data/MARCS                       (15.5k MARCS .mod models)
    av3d     -> storm/data/stagger_avg_mul23_format    (<3D> Stagger averages, MULTI format)
    stagger2 -> storm/data/StaggerGrid2/Lagae/Stagger_grid   (80x80x240 snapshots)

dispatch.x writes run <set>/<name>.nml to <DATADIR>/<set>/<name>/, so the namelist
argument stays a short relative path (SUBAGENTS.md) while the outputs land under runs/.

    python3 make_nml.py bench            # the phase-0 timing set, into bench/
    python3 make_nml.py batch1           # 1D LTE grid over the 4MOST domain, into b1/
    python3 make_nml.py pilot [K]        # paired 1D/<3D>/3D LTE+NLTE pilot, K snapshots, into p1/
"""
import re
import sys
from pathlib import Path

CAMPAIGN = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign')
DATADIR = CAMPAIGN / 'runs'

# FGK line list: VALD atoms (ABO broadening), the dedicated hydrogen file and the
# molecules with bands in the 4MOST windows above ~4000 K.  TiO, VO, ZrO and H2O
# (0.9 GB of the 1 GB master folder) are left out; add TiO for Teff < 4200 K giants.
FGK_LINELISTS = ['vald_atoms_abo.txt', 'Hlinedata',
                 '12CH.txt', '13CH.txt', '12C14N.txt', '13C14N.txt', '12C12C.txt',
                 '12C13C.txt', '12C16O.txt', '13C16O.txt', '16OH.txt', '14NH.txt',
                 'MgH.txt', '28SiH.txt', 'CaH.txt', 'FeH.txt']
MASTER = 'input/linelists/master_linelists'


def group(name, **keys):
    def fmt(v):
        if isinstance(v, bool):
            return '.true.' if v else '.false.'
        if isinstance(v, str):
            return f"'{v}'"
        if isinstance(v, (list, tuple)):          # one value per line: short records
            return ',\n    '.join(fmt(x) for x in v)
        return repr(v)
    body = ' '.join(f'{k}={fmt(v)}' for k, v in keys.items())
    return f'&{name} {body} /'


def linelist_group(linelist, exclude=''):
    keys = {}
    if linelist == 'full':
        keys['linelist_folder'] = MASTER
    elif linelist == 'fgk':
        keys['linelist_file'] = [f'{MASTER}/{f}' for f in FGK_LINELISTS]
    else:                                    # a single (pruned) file
        keys['linelist_file'] = linelist
    if exclude:
        keys['exclude_elements'] = exclude
    return group('linelist_params', **keys)


def namelist(atmos, window='windows/hrs.txt', R=2e5, linelist='fgk', geometry='1d',
             atom=None, feh=None, output=None, rt_output=None, composition=None, extra=()):
    """One run.  geometry '1d' applies the 1d preset (nz=64, exit planes in memory),
    '3d' the 3d preset (30x30x96, tabulated opacities, rotated rays, exit planes on
    disk).  atom: a model atom for NLTE, solved and synthesised in one run (nlte
    preset), with its element excluded from the line list.  output / rt_output: keys
    added to &output_params / &rt_output_params (a group counts only once)."""
    lines = [group('io_params', datadir=str(DATADIR), do_trace=False),
             group('timer_params', sec_per_report=1e8, verbose=2)]
    a = dict(atmos)
    if feh is not None:
        a['FeH'] = feh
    if geometry == '1d':
        a.setdefault('nz', 64)
        lines += [group('atmos_params', **a),
                  group('output_params', rt_to_file=False, **(output or {}))]
    else:
        a.update(nx=30, ny=30, nz=96)
        lines += [group('atmos_params', **a),
                  group('run_params', rotate_atmos=True),
                  group('output_params', rt_to_file=True, **(output or {}))]
    rt_it, rt_out = {}, dict(rt_output or {})
    if geometry == '3d':
        rt_it.update(sample_back_opac=True, sample_linelist=True)
        rt_out.update(sample_back_opac=True, sample_linelist=True)
    exclude = ''
    if atom:
        element = atom['element']
        exclude = element
        lines.append(group('atom_params', atom_file=atom['file']))
        lines.append(group('nlte_params', convlim=1e-3, **atom.get('nlte', {})))
        lines.append(group('opac_table_params', absmet_file='input/absmet_files/absmet'))
        rt_out['maxiter'] = 3
        if atom.get('departures_only'):
            rt_out['angle_quad'] = 'none'
    if rt_it:
        lines.append(group('rt_iteration_params', **rt_it))
    if rt_out:
        lines.append(group('rt_output_params', **rt_out))
    if linelist and not (atom and atom.get('departures_only')):
        lines.append(linelist_group(linelist, exclude))
        lines.append(group('spectrum_params', window_file=window, R_sampling=R))
    if composition:
        lines.append(group('composition_params', **composition))
    lines += list(extra)
    return '\n'.join(lines) + '\n'


MARCS = {
    'sun': 'marcs/p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod',
    'giant': 'marcs/s4750_g+2.5_m1.0_t02_st_z-0.50_a+0.20_c+0.00_n+0.00_o+0.20_r+0.00_s+0.00.mod',
    'mpdwarf': 'marcs/p6000_g+4.0_m0.0_t01_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod',
}
MG = dict(element='Mg', file='input/atoms/atom.mg86e')
FE = dict(element='Fe', file='input/atoms/atom.fe607a', nlte=dict(convmax=1e-3))
SNAP = 'stagger2/t5777/t5777g44m00/atm3d.t5777g44m0005_00076_80x80x240'


def bench():
    m = lambda k: dict(atmos_format='Marcs', atmos_file=MARCS[k])
    runs = {
        # 1D LTE, thread scaling on the HRS windows (64k points at R=2e5)
        'l1d_sun_hrs': namelist(m('sun')),
        'l1d_sun_hrs_full': namelist(m('sun'), linelist='full'),
        'l1d_sun_hrs_R1e5': namelist(m('sun'), R=1e5),
        'l1d_giant_hrs': namelist(m('giant')),
        'l1d_mpdwarf_hrs': namelist(m('mpdwarf')),
        'l1d_sun_lrs': namelist(m('sun'), window='windows/lrs.txt', R=5e4),
        'l1d_sun_hrs_pruned': namelist(m('sun'), linelist='bench/pruned_sun_hrs.txt'),
        # 16 abundance-offset sets on top of the reference, in one process
        'l1d_sun_hrs_sets': namelist(m('sun'), output=dict(set_file=str(CAMPAIGN / 'jobs/sets/sets_x8.txt'))),
        'l1d_sun_hrs_nolinelist_out': namelist(m('sun'), rt_output=dict(save_linelist=False)),
        # <3D>: tau500-averaged Stagger model in MULTI format
        'av3d_sun_hrs': namelist(dict(atmos_file='av3d/atmos.t5777g44m00'), feh=0.0),
        # 1D NLTE, departures and the HRS spectrum in one run
        'n1d_mg_sun_hrs': namelist(m('sun'), atom=MG),
        'n1d_fe_sun_hrs': namelist(m('sun'), atom=FE),
        'n1d_mg_giant_hrs': namelist(m('giant'), atom=MG),
        # 3D LTE on one Stagger-grid snapshot, HRS at R=1e5
        'l3d_sun_hrs_R1e5': namelist(dict(atmos_format='stagger2', atmos_file=SNAP),
                                     geometry='3d', R=1e5, feh=0.0),
        # 3D NLTE Mg departures only (no spectrum)
        'n3d_mg_sun_dep': namelist(dict(atmos_format='stagger2', atmos_file=SNAP),
                                   geometry='3d', feh=0.0,
                                   atom=dict(MG, departures_only=True,
                                             nlte=dict(ng_step=20))),
    }
    # hop_prune.x input for the pruned 1D run: the run's own groups plus &prune_params
    prune = namelist(m('sun')).replace(
        f"{group('timer_params', sec_per_report=1e8, verbose=2)}\n", '')
    prune += group('prune_params', output_file='bench/pruned_sun_hrs.txt', sample_count=8) + '\n'
    out = Path('bench')
    out.mkdir(exist_ok=True)
    for name, text in runs.items():
        (out / f'{name}.nml').write_text(text)
    (out / 'prune_sun_hrs.nml').write_text(prune)
    print(f'wrote {len(runs) + 1} namelists to {out}/')


MARCS_NAME = re.compile(r'(?P<geo>[ps])(?P<teff>\d+)_g(?P<logg>[+-]\d\.\d)_m(?P<mass>[\d.]+)_'
                        r't(?P<vt>\d\d)_st_z(?P<feh>[+-]\d\.\d\d)_a(?P<alpha>[+-]\d\.\d\d)_')


def marcs_readable(path):
    """Whether both depth tables of a MARCS model parse as the Fortran reader expects:
    rows numbered 1..nz and whitespace-separated numbers.  39 domain models (mostly
    hot spherical models at log g <= 2.5) have fused columns and abort the run with
    'model and loop indices differ'."""
    lines = Path(path).read_text().splitlines()
    i = next(j for j, line in enumerate(lines) if 'Number of depth points' in line)
    nz = int(lines[i].split()[0])
    for start in (i + 3, i + 4 + nz):
        for k in range(nz):
            words = lines[start + k].split()
            try:
                if int(words[0]) != k + 1:
                    return False
                [float(w) for w in words]
            except (ValueError, IndexError):
                return False
    return True


def domain_models(folder='marcs'):
    """Standard-composition MARCS models of the 4MOST FGK domain: 4000-6750 K,
    log g 1.0-5.0, [Fe/H] -4.0..+0.5; plane-parallel (vt=1, m=0) at log g >= 3.5 and
    spherical (vt=2, 1 Msun) at log g <= 3.0, so each node has one model."""
    out = []
    for path in sorted(Path(folder).iterdir()):
        m = MARCS_NAME.match(path.name)
        if not m:
            continue
        teff, logg, feh = int(m['teff']), float(m['logg']), float(m['feh'])
        if not (4000 <= teff <= 6750 and 1.0 <= logg <= 5.0 and -4.0 <= feh <= 0.5):
            continue
        plane = m['geo'] == 'p' and m['vt'] == '01' and m['mass'] == '0.0' and logg >= 3.5
        sphere = m['geo'] == 's' and m['vt'] == '02' and m['mass'] == '1.0' and logg <= 3.0
        if (plane or sphere) and marcs_readable(path):
            out.append((path, float(m['alpha'])))
    return out


def batch1():
    """1D LTE over the domain: per model the HRS windows at R=2e5 with the eight-element
    abundance ladder (sets_x8.txt) and line diagnostics, and the LRS range at R=5e4.
    The alpha enhancement of each MARCS model is applied to the magg22 abundances, so
    the synthesis composition follows the atmosphere's."""
    sets = str(CAMPAIGN / 'jobs/sets/sets_x8.txt')
    tasks = []
    for path, alpha in domain_models():
        name = path.name.split('_st_')[0] + '_z' + path.name.split('_z')[1][:5]
        atmos = dict(atmos_format='Marcs', atmos_file=f'marcs/{path.name}')
        comp = dict(alpha_enhancement=alpha) if alpha else None
        for arm, text in (
                ('hrs', namelist(atmos, composition=comp, output=dict(set_file=sets))),
                ('lrs', namelist(atmos, window='windows/lrs.txt', R=5e4, composition=comp,
                                 rt_output=dict(save_linelist=False)))):
            nml = Path('b1') / arm / f'{name}.nml'
            nml.parent.mkdir(parents=True, exist_ok=True)
            nml.write_text(text)
            tasks.append(str(nml))
    # interleave cheap and expensive tasks so array shards finish together
    tasks.sort(key=lambda t: ('hrs' in t, t))
    Path('b1/tasks.txt').write_text('\n'.join(tasks) + '\n')
    print(f'wrote {len(tasks)} namelists and b1/tasks.txt')


# Pilot models in the 1D/3D overlap: Stagger-grid folder, MARCS twin at the nominal
# parameters, <3D> Stagger average, [Fe/H], [alpha/Fe].  The two 4MOST populations
# (turn-off dwarfs near 6000/4.0, red clump/RGB near 4700/2.5) at solar and -2 metallicity,
# plus the Sun and a cooler dwarf.
PILOT = [
    ('t5777/t5777g44m00', 'input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod', 't5777g44m00', 0.0, 0.0),
    ('t60/t60g40m00', 'marcs/p6000_g+4.0_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod', 't60g40m00', 0.0, 0.0),
    ('t60/t60g40m20', 'marcs/p6000_g+4.0_m0.0_t01_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod', 't60g40m20', -2.0, 0.4),
    ('t55/t55g45m00_new', 'marcs/p5500_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod', 't55g45m00', 0.0, 0.0),
    ('t50/t50g25m00', 'marcs/s5000_g+2.5_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod', 't50g25m00', 0.0, 0.0),
    ('t45/t45g25m00_new', 'marcs/s4500_g+2.5_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod', 't45g25m00', 0.0, 0.0),
    ('t45/t45g25m20_new', 'marcs/s4500_g+2.5_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod', 't45g25m20', -2.0, 0.4),
]
# tumult.synthesis.snapshot_grid.halton_order(10): a prefix of k is the k-snapshot pick
HALTON_10 = [7, 2, 9, 5, 4, 8, 0, 6, 3, 1]
FE_PILOT = dict(element='Fe', file='input/atoms/atom.fe607a', nlte=dict(convmax=1e-3))


def pilot(k=2):
    """Paired calculations per pilot model, all on the HRS windows at R=2e5 with the
    same line list and composition: 1D LTE / NLTE (Mg, Fe) on the MARCS twin and on the
    <3D> average, and per snapshot 3D LTE and 3D NLTE (Mg; Fe separately, it needs
    ~80 GB and hours per snapshot).  NLTE runs solve the departures and write the NLTE
    and LTE spectra in one run (nlte preset)."""
    tasks = {'1d': [], '3d_lte': [], '3d_mg': [], '3d_fe': []}

    def put(kind, rel, text):
        nml = Path('p1') / rel
        nml.parent.mkdir(parents=True, exist_ok=True)
        nml.write_text(text)
        tasks[kind].append(str(nml))

    for folder, marcs, av, feh, alpha in PILOT:
        model = Path(folder).name.replace('_new', '')
        comp = dict(alpha_enhancement=alpha) if alpha else None
        one_d = {'marcs': dict(atmos_format='Marcs', atmos_file=marcs),
                 'av3d': dict(atmos_file=f'av3d/atmos.{av}', FeH=feh)}
        for kind, atmos in one_d.items():
            put('1d', f'{model}/{kind}_lte.nml', namelist(atmos, composition=comp))
            put('1d', f'{model}/{kind}_mg.nml', namelist(atmos, atom=MG, composition=comp))
            put('1d', f'{model}/{kind}_fe.nml', namelist(atmos, atom=FE_PILOT, composition=comp))
        snaps = sorted((Path('stagger2') / folder).glob('atm3d.*'))
        for i in sorted(HALTON_10[:k]):
            snap = snaps[i]
            tag = snap.name.split('_')[1]                       # snapshot number
            atmos = dict(atmos_format='stagger2', atmos_file=str(snap))
            put('3d_lte', f'{model}/s{tag}_lte.nml',
                namelist(atmos, geometry='3d', feh=feh, composition=comp))
            put('3d_mg', f'{model}/s{tag}_mg.nml',
                namelist(atmos, geometry='3d', feh=feh, atom=dict(MG, nlte=dict(ng_step=20)),
                         composition=comp))
            put('3d_fe', f'{model}/s{tag}_fe.nml',
                namelist(atmos, geometry='3d', feh=feh,
                         atom=dict(FE_PILOT, nlte=dict(convmax=1e-3, ng_step=20)),
                         composition=comp))
    for kind, names in tasks.items():
        Path(f'p1/tasks_{kind}.txt').write_text('\n'.join(names) + '\n')
    print(', '.join(f'{len(v)} {k}' for k, v in tasks.items()), 'namelists in p1/')


if __name__ == '__main__':
    if sys.argv[1:2] == ['bench']:
        bench()
    elif sys.argv[1:2] == ['batch1']:
        batch1()
    elif sys.argv[1:2] == ['pilot']:
        pilot(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    else:
        sys.exit(__doc__)
