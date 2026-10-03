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
    python3 make_nml.py batch2           # 1D NLTE Mg/Fe map and <3D>-vs-MARCS pairs, into b2/
    python3 make_nml.py batch2_sed       # stage 2: line-list NLTE SEDs on the b2 departures
    python3 make_nml.py pilot_sed [N]    # stage 2 for the pilot: 1D, 3D and mapped-<3D> 3D SEDs (N output scattering iterations)
    python3 make_nml.py vmic_av3d        # <3D> LTE with the MARCS twin's constant vmic, into b3/
    python3 make_nml.py lte3d LABELS [LIST]  # 3D LTE, first snapshot of each Stagger model near 4MOST stars, into x3d/
    python3 make_nml.py avg3d [MODEL ...]     # own <3D> averages with density: snapshot runs, 1D LTE/Mg/Fe, into avg3d/ and b4/
    python3 make_nml.py avg3d_sed MODEL      # stage-2 SED on a model's b4 departures, once they exist
    python3 make_nml.py vmic_grid            # 1D MARCS and <3D> LTE on a vmic grid for every 3D LTE model, into b5/
    python3 make_nml.py sed_validate         # cost/accuracy test of the 3D SED output-pass settings, into val/
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
    preset), with its element excluded from the line list; atom['nlte'] adds to
    &nlte_params and atom['maxiter'] sets the iteration limit.  output / rt_output: keys
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
        lines.append(group('nlte_params', **{'convlim': 1e-3, **atom.get('nlte', {})}))
        lines.append(group('opac_table_params', absmet_file='input/absmet_files/absmet'))
        rt_out['maxiter'] = 3
        if 'maxiter' in atom:
            rt_it['maxiter'] = atom['maxiter']
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
# One stopping rule per element for every approximation (1D, <3D>, 3D), so the paired
# spectra share it.  convlim is the global criterion; convmax = convlim/10 only freezes
# patches that are well inside it.  With convmax = convlim the patches end the run one
# iteration before the global change confirms convlim, with twice the flux error of
# convlim alone.  A run that ends on convmax quotes the previous iteration's global change
# and can warn of non-convergence although every patch changed by less than convmax.
# In 1D against convlim=1e-5 references, convlim=1e-3 leaves up to 7e-4 in the
# normalised flux, 1e-4 at most 4e-5 (reports/2026-10-03_phase1_launch.md).
# Mg is cheap and takes 1e-4; Fe takes 1e-3, which the weekend's 3D budget allows.  The
# 3D runs converge more slowly per iteration (Mg: ~0.92 against ~0.8 in 1D), hence
# Ng every 20 iterations (as in the phase-0 3D Mg run) and a higher iteration limit.
PILOT_MG = dict(element='Mg', file='input/atoms/atom.mg86e', nlte=dict(convlim=1e-4, convmax=1e-5))
PILOT_FE = dict(element='Fe', file='input/atoms/atom.fe607a', nlte=dict(convlim=1e-3, convmax=1e-4))


def in_3d(atom):
    return dict(atom, nlte=dict(atom['nlte'], ng_step=20), maxiter=150)


def pilot(k=2):
    """Paired calculations per pilot model, all on the HRS windows at R=2e5 with the
    same line list and composition: 1D LTE / NLTE (Mg, Fe) on the MARCS twin and on the
    <3D> average, and per snapshot 3D LTE and 3D NLTE (Mg; Fe separately, it needs
    ~80 GB and hours per snapshot).  NLTE runs solve the departures and write the NLTE
    and LTE spectra in one run (nlte preset); that LTE spectrum shares the model atom's
    lines, so Delta_NLTE, Delta_3D and Delta_interaction of an element come from its NLTE
    runs alone.  The plain LTE runs use the VALD lines of every element.

    The 3D task lists are ordered by snapshot rank (the first Halton snapshot of every
    model before the second), so the first array elements cover all models."""
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
            put('1d', f'{model}/{kind}_mg.nml', namelist(atmos, atom=PILOT_MG, composition=comp))
            put('1d', f'{model}/{kind}_fe.nml', namelist(atmos, atom=PILOT_FE, composition=comp))
    for i in HALTON_10[:k]:
        for folder, marcs, av, feh, alpha in PILOT:
            model = Path(folder).name.replace('_new', '')
            comp = dict(alpha_enhancement=alpha) if alpha else None
            snap = sorted((Path('stagger2') / folder).glob('atm3d.*'))[i]
            tag = snap.name.split('_')[1]                       # snapshot number
            atmos = dict(atmos_format='stagger2', atmos_file=str(snap))
            put('3d_lte', f'{model}/s{tag}_lte.nml',
                namelist(atmos, geometry='3d', feh=feh, composition=comp))
            put('3d_mg', f'{model}/s{tag}_mg.nml',
                namelist(atmos, geometry='3d', feh=feh, atom=in_3d(PILOT_MG), composition=comp))
            put('3d_fe', f'{model}/s{tag}_fe.nml',
                namelist(atmos, geometry='3d', feh=feh, atom=in_3d(PILOT_FE), composition=comp))
    for kind, names in tasks.items():
        Path(f'p1/tasks_{kind}.txt').write_text('\n'.join(names) + '\n')
    print(', '.join(f'{len(v)} {k}' for k, v in tasks.items()), 'namelists in p1/')


AV3D_NAME = re.compile(r't(?P<teff>\d+)g(?P<logg>\d\d)(?P<sign>[mp])(?P<feh>\d\d)$')


def av3d_labels(name):
    """Nominal Teff, log g and [Fe/H] of a Stagger <3D> average from its name
    (t45g25m05: 4500 K, 2.5, -0.5; t5777g44p05: 5777 K, 4.4, +0.5), and the alpha
    enhancement of the Stagger grid: +0.4 at [Fe/H] <= -1, otherwise 0 (Magic et al.
    2013)."""
    m = AV3D_NAME.match(name)
    teff = int(m['teff']) * (100 if len(m['teff']) == 2 else 1)
    feh = int(m['feh']) / 10 * (-1 if m['sign'] == 'm' else 1) + 0.0     # no -0.0
    return teff, int(m['logg']) / 10, feh, 0.4 if feh <= -1 else 0.0


def marcs_twin(teff, logg, feh, folder='marcs'):
    """The standard-composition MARCS model at the node of a Stagger average, under
    the batch-1 geometry rule (plane-parallel at log g >= 3.5, spherical 1 Msun at
    log g <= 3.0).  The Sun takes the dedicated p5777_g+4.4 model at solar
    metallicity and p5750_g+4.5 otherwise, the nearest MARCS node."""
    if (teff, logg) == (5777, 4.4):
        if feh == 0:
            return Path(PILOT[0][1])
        teff, logg = 5750, 4.5
    geo, mass, vt = ('p', '0.0', '01') if logg >= 3.5 else ('s', '1.0', '02')
    hits = sorted(Path(folder).glob(f'{geo}{teff}_g{logg:+.1f}_m{mass}_t{vt}_st_z{feh:+.2f}_*.mod'))
    return hits[0] if hits and marcs_readable(hits[0]) else None


def run_name(path):
    """batch-1 run name of a MARCS model: p5750_g+4.5_m0.0_t01_z+0.00."""
    return path.name.split('_st_')[0] + '_z' + path.name.split('_z')[1][:5]


def batch2():
    """1D NLTE Mg and Fe over the batch-1 MARCS models (A) and LTE, NLTE Mg and NLTE Fe
    on every Stagger <3D> average with its MARCS twin (B), all on the HRS windows at
    R=2e5 with the pilot's convergence rules.  Each NLTE run writes the NLTE and LTE
    spectra of its element and keeps the departure coefficients.  linelist_out is not
    written (it would be 80 % of the output).

    A twin is taken from A (and batch 1 for LTE) when the model and the synthesis
    composition are the same; the Sun, 7000 K and the [Fe/H] = -0.5 averages (Stagger
    alpha 0, MARCS +0.2) get their own twin runs with the <3D> composition.
    b2/pairs.txt lists every <3D> run folder with its labels and twin folders."""
    quiet = dict(save_linelist=False)
    tasks = {'a': [], 'b': []}

    def put(kind, rel, text):
        nml = Path('b2') / rel
        nml.parent.mkdir(parents=True, exist_ok=True)
        nml.write_text(text)
        tasks[kind].append(str(nml))
        return str(Path('b2') / Path(rel).with_suffix(''))

    in_a = {}
    for path, alpha in domain_models():
        name = run_name(path)
        atmos = dict(atmos_format='Marcs', atmos_file=f'marcs/{path.name}')
        comp = dict(alpha_enhancement=alpha) if alpha else None
        for el, atom in (('mg', PILOT_MG), ('fe', PILOT_FE)):
            put('a', f'n1d/{name}_{el}.nml', namelist(atmos, atom=atom, composition=comp,
                                                      rt_output=quiet))
        in_a[path.name] = alpha
    pairs = ['# av3d teff logg feh alpha  av3d_run  twin_lte twin_mg twin_fe  twin_model']
    for av in sorted(p.name[6:] for p in Path('av3d').glob('atmos.t*')):
        teff, logg, feh, alpha = av3d_labels(av)
        comp = dict(alpha_enhancement=alpha) if alpha else None
        atmos = dict(atmos_file=f'av3d/atmos.{av}', FeH=feh)
        put('b', f'av3d/{av}/lte.nml', namelist(atmos, composition=comp, rt_output=quiet))
        for el, atom in (('mg', PILOT_MG), ('fe', PILOT_FE)):
            put('b', f'av3d/{av}/{el}.nml', namelist(atmos, atom=atom, composition=comp,
                                                     rt_output=quiet))
        twin = marcs_twin(teff, logg, feh)
        if twin is None:
            runs = ['-'] * 3
        elif in_a.get(twin.name) == alpha:
            runs = [f'b1/hrs/{run_name(twin)}', f'b2/n1d/{run_name(twin)}_mg',
                    f'b2/n1d/{run_name(twin)}_fe']
        else:
            m = dict(atmos_format='Marcs', atmos_file=str(twin))
            runs = [put('b', f'twin/{av}/lte.nml', namelist(m, composition=comp, rt_output=quiet))]
            runs += [put('b', f'twin/{av}/{el}.nml', namelist(m, atom=atom, composition=comp,
                                                              rt_output=quiet))
                     for el, atom in (('mg', PILOT_MG), ('fe', PILOT_FE))]
        pairs.append(f'{av} {teff} {logg} {feh} {alpha}  b2/av3d/{av}  {" ".join(runs)}  '
                     f'{twin.name if twin else "-"}')
    # the expensive Fe runs first, so the array's tail holds the short Mg and LTE runs
    order = lambda t: (not t.endswith('fe.nml'), t)
    for kind, names in tasks.items():
        Path(f'b2/tasks_{kind}.txt').write_text('\n'.join(sorted(names, key=order)) + '\n')
    Path('b2/pairs.txt').write_text('\n'.join(pairs) + '\n')
    print(', '.join(f'{len(v)} {k}' for k, v in tasks.items()), 'namelists in b2/,',
          sum(' - ' in p for p in pairs), '<3D> averages without a twin')


def combine_departures(folder, mg_run, fe_run):
    """A read_departures source holding the Mg and Fe departures of two runs on the
    same atmosphere: params.nml, atmos_mesh.txt and save/* of the Fe run, plus
    save/pop_Mg_* of the Mg run, as symlinks."""
    save = folder / 'save'
    save.mkdir(parents=True, exist_ok=True)
    links = [(fe_run / f, folder / f) for f in ('params.nml', 'atmos_mesh.txt')]
    links += [(f, save / f.name) for f in (fe_run / 'save').iterdir()]
    links += [(f, save / f.name) for f in (mg_run / 'save').glob('pop_Mg_*')]
    for src, dst in links:
        if not dst.is_symlink():
            dst.symlink_to(src)


def sed_namelist(dep_folder, composition=None, window='windows/hrs.txt', R=2e5):
    """Stage 2 of the NLTE synthesis: the spectrum from the full FGK line list on the
    precomputed Mg and Fe departures.  The atoms keep their levels and continua but not
    their lines (disable_bb), and their species stay in the line list
    (exclude_from_line_list=F), so every Mg and Fe line of the list takes the departure
    coefficients of the atom levels nearest in energy.  The run's LTE spectrum
    (crop(LTE=True)) is the matching LTE SED: same line list, mesh and wavelengths."""
    lines = [group('io_params', datadir=str(DATADIR), do_trace=False),
             group('timer_params', sec_per_report=1e8, verbose=2),
             group('atmos_params', atmos_format='save', atmos_file=str(dep_folder),
                   read_departures=True),
             group('output_params', rt_to_file=False),
             group('atom_params', atom_file=[PILOT_MG['file'], PILOT_FE['file']],
                   disable_bb=True, exclude_from_line_list=False),
             group('opac_table_params', absmet_file='input/absmet_files/absmet'),
             group('rt_iteration_params', maxiter=0),
             group('rt_output_params', maxiter=3, save_linelist=False),
             linelist_group('fgk'),
             group('spectrum_params', window_file=window, R_sampling=R)]
    if composition:
        lines.append(group('composition_params', **composition))
    return '\n'.join(lines) + '\n'


def sed_namelist_3d(dep_folder, atoms, feh, composition=None):
    """Stage 2 on a 3D snapshot: as sed_namelist, with the 3d preset's tabulated
    opacities and rotated rays, for the given model atoms."""
    text = sed_namelist(dep_folder, composition)
    text = text.replace(group('output_params', rt_to_file=False),
                        group('run_params', rotate_atmos=True) + '\n' +
                        group('output_params', rt_to_file=True))
    text = text.replace(group('rt_iteration_params', maxiter=0),
                        group('rt_iteration_params', maxiter=0, sample_back_opac=True,
                              sample_linelist=True))
    text = text.replace(group('rt_output_params', maxiter=3, save_linelist=False),
                        group('rt_output_params', maxiter=3, save_linelist=False,
                              sample_back_opac=True, sample_linelist=True))
    return text.replace(group('atom_params', atom_file=[PILOT_MG['file'], PILOT_FE['file']],
                              disable_bb=True, exclude_from_line_list=False),
                        group('atom_params', atom_file=atoms, disable_bb=True,
                              exclude_from_line_list=False))


def solved(run, element):
    """Whether an NLTE run has written the departures of all its depth patches.  Each
    patch saves them (one row of pop_<el>_meta.txt) when it stops iterating, before the
    output pass (about 3 h of a 3D run), so a run that its time limit ends during the
    output pass still serves stage 2."""
    rows = lambda f: (sum(1 for l in f.read_text().splitlines() if not l.startswith('*'))
                      if f.is_file() else -1)
    save = run / 'save'
    return rows(save / f'pop_{element}_meta.txt') == rows(save / 'atmos_meta.txt')


def pilot_sed(rt_maxiter=3):
    """Stage-2 SEDs of the pilot, Mg departures (Mg + Fe once a 3D Fe run has finished), all from the
    full line list: 1D on the MARCS twin and the <3D> average (b2-style combined Mg+Fe
    departures), and per finished 3D Mg snapshot (a) on its own 3D departures and (b) on
    the <3D> departures mapped to it in log tau_500 (map_departures.py, run by the job
    before the SED).  The mapped departures come from the own <3D> averages with density
    (b4/av3d/<model>, make_nml.py avg3d) into runs/p1sed/map4/<model>_s<snap>; the
    earlier maps from the Stagger files (runs/p1sed/map/) stay as they are.

    p1sed/tasks_3d.txt is in priority order: true-3D-departure SEDs, then the mapped
    ones of each model's first snapshot, then those of the second (snapshot rank as
    in pilot(): HALTON_10 over the sorted snapshots).  rt_maxiter sets the 3D SEDs'
    &rt_output_params maxiter (scattering iterations of the output pass).  Writes
    p1sed/tasks_1d.txt, p1sed/tasks_3d.txt and p1sed/map.txt (1D run, 3D template, out)."""
    t1, t3, maps = [], [], []

    def sed3d(dep, atoms, feh, comp):
        text = sed_namelist_3d(dep, atoms, feh, comp)
        return text.replace('maxiter=3 save_linelist', f'maxiter={rt_maxiter} save_linelist')

    for folder, marcs, av, feh, alpha in PILOT:
        model = Path(folder).name.replace('_new', '')
        comp = dict(alpha_enhancement=alpha) if alpha else None
        snaps = [p.name.split('_')[1] for p in sorted((Path('stagger2') / folder).glob('atm3d.*'))]
        rank = {f's{tag}': HALTON_10.index(i) if i in HALTON_10 else 99 for i, tag in enumerate(snaps)}
        for kind in ('marcs', 'av3d'):
            mg, fe = DATADIR / f'p1/{model}/{kind}_mg', DATADIR / f'p1/{model}/{kind}_fe'
            dep = DATADIR / 'p1sed/dep' / f'{model}_{kind}'
            combine_departures(dep, mg, fe)
            nml = Path('p1sed') / model / f'{kind}.nml'
            nml.parent.mkdir(parents=True, exist_ok=True)
            nml.write_text(sed_namelist(dep, comp))
            t1.append(str(nml))
        av_new = DATADIR / f'b4/av3d/{model}'
        for run in sorted((DATADIR / f'p1/{model}').glob('s*_mg')):
            if not solved(run, 'Mg'):
                continue                                   # still iterating
            snap = run.name.split('_')[0]
            nml = Path('p1sed') / model / f'{snap}_3dmg.nml'
            nml.write_text(sed3d(run, [PILOT_MG['file']], feh, comp))
            t3.append((0, rank[snap], str(nml)))
            if solved(av_new / 'mg', 'Mg'):
                out = DATADIR / 'p1sed/map4' / f'{model}_{snap}'
                maps.append(f'{av_new}/mg {DATADIR}/p1/{model}/{snap}_lte {out} Mg')
                nml = Path('p1sed') / model / f'{snap}_map4mg.nml'
                nml.write_text(sed3d(out, [PILOT_MG['file']], feh, comp))
                t3.append((1 + rank[snap], 1, str(nml)))
            fe = DATADIR / f'p1/{model}/{snap}_fe'
            if solved(fe, 'Fe'):
                dep = DATADIR / 'p1sed/dep' / f'{model}_{snap}'
                combine_departures(dep, run, fe)
                nml = Path('p1sed') / model / f'{snap}_3dmgfe.nml'
                nml.write_text(sed3d(dep, [PILOT_MG['file'], PILOT_FE['file']], feh, comp))
                t3.append((0, rank[snap], str(nml)))
                if solved(av_new / 'mg', 'Mg') and solved(av_new / 'fe', 'Fe'):
                    out = DATADIR / 'p1sed/map4' / f'{model}_{snap}_mgfe'
                    for el in ('mg', 'fe'):                        # both maps into one folder
                        maps.append(f'{av_new}/{el} {DATADIR}/p1/{model}/{snap}_lte '
                                    f'{out} {el.capitalize()}')
                    nml = Path('p1sed') / model / f'{snap}_map4mgfe.nml'
                    nml.write_text(sed3d(out, [PILOT_MG['file'], PILOT_FE['file']], feh, comp))
                    t3.append((1 + rank[snap], 0, str(nml)))
    Path('p1sed/tasks_1d.txt').write_text('\n'.join(t1) + '\n')
    Path('p1sed/tasks_3d.txt').write_text('\n'.join(n for *_, n in sorted(t3)) + '\n')
    Path('p1sed/map.txt').write_text('\n'.join(maps) + '\n')
    print(f'{len(t1)} 1D and {len(t3)} 3D stage-2 namelists in p1sed/')


def batch2_sed():
    """Stage-2 NLTE SEDs of every batch-2 departure pair (A: MARCS grid, B: <3D>
    averages and their own twins), needing the f63e2fc95 build (disable_bb for classic
    atoms, several precomputed atoms).  Departure folders: runs/b2/dep/<run>."""
    tasks = []

    def one(rel, mg, fe, comp):
        dep = DATADIR / 'b2/dep' / rel
        if not ((DATADIR / mg / 'save').is_dir() and (DATADIR / fe / 'save').is_dir()):
            print(f'missing departures for {rel}')
            return
        combine_departures(dep, DATADIR / mg, DATADIR / fe)
        nml = Path('b2/sed') / f'{rel}.nml'
        nml.parent.mkdir(parents=True, exist_ok=True)
        nml.write_text(sed_namelist(dep, comp))
        tasks.append(str(nml))

    for path, alpha in domain_models():
        name = run_name(path)
        one(name, f'b2/n1d/{name}_mg', f'b2/n1d/{name}_fe',
            dict(alpha_enhancement=alpha) if alpha else None)
    for row in Path('b2/pairs.txt').read_text().splitlines()[1:]:
        av, _, _, _, alpha, av_run, _, tw_mg, tw_fe = row.split()[:9]
        comp = dict(alpha_enhancement=float(alpha)) if float(alpha) else None
        one(f'av3d_{av}', f'{av_run}/mg', f'{av_run}/fe', comp)
        if tw_mg.startswith('b2/twin/'):
            one(f'twin_{av}', tw_mg, tw_fe, comp)
    Path('b2/tasks_sed.txt').write_text('\n'.join(tasks) + '\n')
    print(f'{len(tasks)} stage-2 namelists in b2/sed/')


def vmic_av3d():
    """The batch-2 <3D> LTE runs again with the constant microturbulence of their MARCS
    twin (1 km/s at log g >= 3.5, 2 km/s below, as the twin's t01/t02), so that a
    <3D>-MARCS difference holds the temperature and density structure only.  The <3D>
    files carry their own depth-dependent vmic (about 1.4-2.5 km/s in giants).  Into
    b3/av3d_vmic/<av>/lte.nml."""
    tasks = []
    for row in Path('b2/pairs.txt').read_text().splitlines()[1:]:
        av, teff, logg, feh, alpha = row.split()[:5]
        comp = dict(alpha_enhancement=float(alpha)) if float(alpha) else None
        vmic = 1.0 if float(logg) >= 3.5 else 2.0
        atmos = dict(atmos_file=f'av3d/atmos.{av}', FeH=float(feh), vmic=vmic)
        nml = Path('b3/av3d_vmic') / av / 'lte.nml'
        nml.parent.mkdir(parents=True, exist_ok=True)
        nml.write_text(namelist(atmos, composition=comp, rt_output=dict(save_linelist=False)))
        tasks.append(str(nml))
    Path('b3/tasks_vmic.txt').write_text('\n'.join(tasks) + '\n')
    print(f'{len(tasks)} <3D> LTE namelists with the twin vmic in b3/av3d_vmic/')


def stagger_models(root='stagger2'):
    """{<3D> name: Stagger-grid folder} of every 3D model, the '_new' run where a model
    has two (t40g25m00)."""
    out = {}
    for d in sorted(Path(root).glob('t*/t*')):
        name = d.name.replace('_new', '')
        if AV3D_NAME.match(name) and (name not in out or d.name.endswith('_new')):
            out[name] = d
    return out


def star_labels(csv_path):
    """(Teff, log g, [Fe/H]) of the 4MOST HRS stars.  From s4_labels.py's table: the
    Gaia XP XGBoost labels (Andrae et al. 2023; [M/H] taken as [Fe/H]) of every star
    that has them.  From place_stars.py's labels.csv (the phase-2 crude 1D LTE fit): one
    point per star, its highest-S/N spectrum, log g held at Gaia where possible, without
    the 0201 stars pinned at the grid edge [Fe/H] = -4."""
    import csv
    rows = list(csv.DictReader(open(csv_path)))
    if rows and 'label_src' in rows[0]:
        return [(float(r['teff']), float(r['logg']), float(r['mh'])) for r in rows if r['label_src']]
    best = {}
    for row in rows:
        star = row['OBJECT'] or row['specid']
        if star not in best or float(row['snr']) > float(best[star]['snr']):
            best[star] = row
    out = []
    for row in best.values():
        held = row['teff_c'] != ''
        t, g, m = (float(row[k + ('_c' if held else '')]) for k in ('teff', 'logg', 'feh'))
        if row['set'] in ('0201', '201') and m <= -3.99:
            continue
        out.append((t, g, m))
    return out


def lte3d(labels_csv, out='x3d/tasks.txt'):
    """3D LTE on the Halton-first snapshot of every Stagger model in the 4MOST domain
    with placed stars within one grid step (250 K, 0.5 dex in log g and [Fe/H]), as the
    pilot's 3D LTE runs (3d preset, HRS at R=2e5, FGK line list, the <3D> composition).
    The pilot models are skipped.  x3d/tasks.txt lists them by star count, most first;
    x3d/counts.txt gives the counts.  The <3D> and MARCS LTE partners are batch 2's
    (b2/pairs.txt)."""
    stars = star_labels(labels_csv)
    pilot = {Path(f).name.replace('_new', '') for f, *_ in PILOT}
    ranked = []
    for name, folder in stagger_models().items():
        teff, logg, feh, alpha = av3d_labels(name)
        if not (4000 <= teff <= 6750 and 1.0 <= logg <= 5.0 and -4.0 <= feh <= 0.5):
            continue
        n = sum(abs(t - teff) <= 250 and abs(g - logg) <= 0.5 and abs(m - feh) <= 0.5
                for t, g, m in stars)
        if n and name not in pilot:
            ranked.append((n, name, folder, feh, alpha))
    ranked.sort(key=lambda r: (-r[0], r[1]))
    tasks, counts = [], []
    for n, name, folder, feh, alpha in ranked:
        snap = sorted(folder.glob('atm3d.*'))[HALTON_10[0]]
        atmos = dict(atmos_format='stagger2', atmos_file=str(snap))
        comp = dict(alpha_enhancement=alpha) if alpha else None
        nml = Path('x3d') / name / f's{snap.name.split("_")[1]}_lte.nml'
        nml.parent.mkdir(parents=True, exist_ok=True)
        nml.write_text(namelist(atmos, geometry='3d', feh=feh, composition=comp))
        tasks.append(str(nml))
        counts.append(f'{name} {n} {nml}')
    Path(out).write_text('\n'.join(tasks) + '\n')
    Path(out).with_name(Path(out).stem.replace('tasks', 'counts') + '.txt').write_text('\n'.join(counts) + '\n')
    print(f'{len(tasks)} 3D LTE namelists in x3d/ ({len(stars)} stars)')


# Snapshots pooled into one <3D> average: evenly spaced over a model's sequence, so the
# average spans its time coverage rather than one stretch of it.
AVG_SNAPSHOTS = 10


def avg3d_models():
    """Models that get their own <3D> average: the pilot, every [Fe/H] <= -3 model of
    batch 2 (b2/pairs.txt) and every model of the 3D LTE extension list."""
    names = [Path(f).name.replace('_new', '') for f, *_ in PILOT]
    for row in Path('b2/pairs.txt').read_text().splitlines()[1:]:
        if float(row.split()[3]) <= -3:
            names.append(row.split()[0])
    tasks = Path('x3d/tasks_s4.txt')
    if tasks.exists():
        names += [Path(t).parent.name for t in tasks.read_text().split()]
    return list(dict.fromkeys(names))


def avg3d(models=()):
    """Own <3D> averages, with density, for the Stagger models of avg3d_models() (or
    the ones named).  Per model:

    - avg3d/<model>/s<snap>.nml: atmosphere-only runs (a 0.5 A window, no line list) on
      AVG_SNAPSHOTS snapshots, through the 3d preset's 30x30x96 resampling and the
      model's composition; each saves T, n_e, rho, velocities and chi_500.
    - avg3d.py pools their columns into runs/../work/avg3d/<model>.atmos (log tau_500
      surfaces; text format).
    - b4/av3d/<model>/{lte,mg,fe}.nml: the batch-2 <3D> runs on that file (Mg and Fe
      for the pilot and [Fe/H] <= -3 only), which a run
      reads with use_ne=F (rho given, n_e from the EOS) and amr='T' (the default
      pressure-based depth sampling needs n_e).

    Writes avg3d/models.txt (model, [Fe/H], the snapshot namelists) for avg3d.sh."""
    rows = []
    folders = stagger_models()
    pilot = {Path(f).name.replace('_new', '') for f, *_ in PILOT}
    missing = []
    for name in (models or avg3d_models()):
        if name not in folders:
            missing.append(name)                    # a <3D> file without snapshots
            continue
        teff, logg, feh, alpha = av3d_labels(name)
        comp = dict(alpha_enhancement=alpha) if alpha else None
        snaps = sorted(folders[name].glob('atm3d.*'))
        pick = sorted({round(i * (len(snaps) - 1) / (AVG_SNAPSHOTS - 1))
                       for i in range(AVG_SNAPSHOTS)}) if len(snaps) > AVG_SNAPSHOTS else range(len(snaps))
        nmls = []
        for i in pick:
            snap = snaps[i]
            nml = Path('avg3d') / name / f's{snap.name.split("_")[1]}.nml'
            nml.parent.mkdir(parents=True, exist_ok=True)
            atmos = dict(atmos_format='stagger2', atmos_file=str(snap), FeH=feh, nx=30, ny=30, nz=96)
            lines = [group('io_params', datadir=str(DATADIR), do_trace=False),
                     group('timer_params', sec_per_report=1e8, verbose=2),
                     group('atmos_params', **atmos),
                     group('spectrum_params', aa_blue=5000.0, aa_red=5000.5, daa=0.1)]
            if comp:
                lines.append(group('composition_params', **comp))
            nml.write_text('\n'.join(lines) + '\n')
            nmls.append(str(nml))
        atmos = dict(atmos_format='text', atmos_file=f'avg3d/{name}.atmos', FeH=feh, use_ne=False, amr='T')
        quiet = dict(save_linelist=False)
        out = Path('b4/av3d') / name
        out.mkdir(parents=True, exist_ok=True)
        (out / 'lte.nml').write_text(namelist(atmos, composition=comp, rt_output=quiet))
        if name in pilot or feh <= -3:            # the 3D LTE extension needs LTE only
            for el, atom in (('mg', PILOT_MG), ('fe', PILOT_FE)):
                (out / f'{el}.nml').write_text(namelist(atmos, atom=atom, composition=comp, rt_output=quiet))
        rows.append(f'{name} {feh} ' + ' '.join(nmls))
    listing = Path('avg3d/models.txt')
    old = {r.split()[0]: r for r in listing.read_text().splitlines()} if listing.exists() else {}
    old.update({r.split()[0]: r for r in rows})
    listing.write_text('\n'.join(old.values()) + '\n')
    print(f'{len(rows)} models in avg3d/models.txt; no snapshots for', ' '.join(missing) or 'none')


def avg3d_sed(name):
    """Stage-2 SED (batch2_sed's settings) on the Mg and Fe departures of the b4 run of
    one own average: b4/sed/av3d_<model>.nml, departures linked in runs/b4/dep/."""
    teff, logg, feh, alpha = av3d_labels(name)
    comp = dict(alpha_enhancement=alpha) if alpha else None
    dep = DATADIR / 'b4/dep' / f'av3d_{name}'
    combine_departures(dep, DATADIR / f'b4/av3d/{name}/mg', DATADIR / f'b4/av3d/{name}/fe')
    nml = Path('b4/sed') / f'av3d_{name}.nml'
    nml.parent.mkdir(parents=True, exist_ok=True)
    nml.write_text(sed_namelist(dep, comp))
    print(nml)


def sed_validate():
    """val/: the cost/accuracy test of the 3D stage-2 SED settings.  The Sun's first
    3D Mg snapshot on its own departures, and the metal-poor giant t45g25m20 on the own
    <3D> Mg departures mapped to a snapshot, each with 3, 1 and 0 scattering iterations
    in the output pass (&rt_output_params maxiter); the reference is maxiter=3, the
    pilot setting.  Writes val/tasks.txt (in run order) and val/map.txt."""
    cases = [('sun', DATADIR / 'p1/t5777g44m00/s00459_mg', 0.0, None, None),
             ('mpg', DATADIR / 'val/map_t45g25m20_s00233', -2.0, dict(alpha_enhancement=0.4),
              f'{DATADIR}/b4/av3d/t45g25m20/mg {DATADIR}/p1/t45g25m20/s00233_lte '
              f'{DATADIR}/val/map_t45g25m20_s00233 Mg')]
    tasks, maps = [], []
    for tag, dep, feh, comp, mapline in cases:
        if mapline:
            maps.append(mapline)
        for it in (3, 1, 0):
            text = sed_namelist_3d(dep, [PILOT_MG['file']], feh, comp)
            nml = Path('val') / f'{tag}_it{it}.nml'
            nml.parent.mkdir(parents=True, exist_ok=True)
            nml.write_text(text.replace('maxiter=3 save_linelist', f'maxiter={it} save_linelist'))
            tasks.append(str(nml))
    order = ['val/sun_it3.nml', 'val/sun_it1.nml', 'val/mpg_it3.nml', 'val/mpg_it1.nml',
             'val/sun_it0.nml', 'val/mpg_it0.nml']
    Path('val/tasks.txt').write_text('\n'.join(order) + '\n')
    Path('val/map.txt').write_text('\n'.join(maps) + '\n')
    print(len(tasks), 'namelists in val/')


VMIC_GRID = (0.5, 1.0, 1.5, 2.0, 2.5)


def vmic_grid():
    """1D LTE partners of every model with a finished 3D LTE snapshot (the pilot's
    p1/<model>/s*_lte and the extension's x3d/<model>/s*_lte), on a constant-vmic grid
    VMIC_GRID: the MARCS twin of b2/pairs.txt and the <3D> average, the own one with
    density where make_nml.py avg3d has written it (avg3d/<model>.atmos), else the
    Stagger file.  Same line list, windows and composition as the 3D runs; into
    b5/vmic/<model>/{marcs,av3d}_v<vmic>.nml.  Rerun as more 3D LTE runs finish; existing
    namelists are rewritten, finished runs are skipped by the queue."""
    pairs = {r.split()[0]: r.split() for r in Path('b2/pairs.txt').read_text().splitlines()[1:]}
    logs = CAMPAIGN / 'logs'
    models = set()
    for top in ('p1', 'x3d'):
        for log in logs.glob(f'{top}_t*_s*_lte.log'):
            if 'OpenMULTI3D terminated' in log.read_text(errors='replace'):
                models.add(log.name[len(top) + 1:].rsplit('_s', 1)[0])
    tasks = []
    for name in sorted(models):
        if name not in pairs:
            continue
        av, teff, logg, feh, alpha, *_ , twin = pairs[name]
        feh, alpha = float(feh), float(alpha)
        comp = dict(alpha_enhancement=alpha) if alpha else None
        own = Path(f'avg3d/{name}.atmos')
        one_d = {'av3d': (dict(atmos_format='text', atmos_file=str(own), FeH=feh, use_ne=False, amr='T')
                          if own.exists() else dict(atmos_file=f'av3d/atmos.{name}', FeH=feh))}
        if twin != '-':
            one_d['marcs'] = dict(atmos_format='Marcs', atmos_file=f'marcs/{twin}')
        for kind, atmos in one_d.items():
            for v in VMIC_GRID:
                nml = Path('b5/vmic') / name / f'{kind}_v{v:.1f}.nml'
                nml.parent.mkdir(parents=True, exist_ok=True)
                nml.write_text(namelist(dict(atmos, vmic=v), composition=comp,
                                        rt_output=dict(save_linelist=False)))
                tasks.append(str(nml))
    Path('b5/tasks_vmic.txt').write_text('\n'.join(tasks) + '\n')
    print(f'{len(tasks)} namelists for {len(models)} models in b5/vmic/')


if __name__ == '__main__':
    if sys.argv[1:2] == ['vmic_grid']:
        vmic_grid()
    elif sys.argv[1:2] == ['sed_validate']:
        sed_validate()
    elif sys.argv[1:2] == ['avg3d']:
        avg3d(sys.argv[2:])
    elif sys.argv[1:2] == ['avg3d_sed']:
        avg3d_sed(sys.argv[2])
    elif sys.argv[1:2] == ['lte3d']:
        lte3d(*sys.argv[2:4])
    elif sys.argv[1:2] == ['vmic_av3d']:
        vmic_av3d()
    elif sys.argv[1:2] == ['bench']:
        bench()
    elif sys.argv[1:2] == ['batch1']:
        batch1()
    elif sys.argv[1:2] == ['batch2']:
        batch2()
    elif sys.argv[1:2] == ['batch2_sed']:
        batch2_sed()
    elif sys.argv[1:2] == ['pilot_sed']:
        pilot_sed(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
    elif sys.argv[1:2] == ['pilot']:
        pilot(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    else:
        sys.exit(__doc__)
