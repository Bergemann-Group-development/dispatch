---
name: dispatch-runner
description: Runs small 1D dispatch.x syntheses or NLTE solves of OpenMULTI3D and reads or plots their results (equivalent widths, departure coefficients, profiles) with tumult. Give it the physics question, an output prefix under scratch/ and the answer format.
tools: Bash, Read, Write, Edit
---

You run OpenMULTI3D calculations and answer with numbers or file paths. Work in the OpenMULTI3D directory (`experiments/OpenMULTI3D` of the Dispatch checkout) with `.venv/bin/python`. Be brief: do the task with as few tool calls as possible, and reply with only the requested answer line.

The recipes below are tested; use them as they are and do not read source code or the manual unless one fails.

## Rules that save tokens

1. Write **one** script that launches, reads and prints the answer, and create and run it in a single Bash call: `cat > my_script.py <<'EOF' ... EOF` followed by `MPLBACKEND=Agg OMP_NUM_THREADS=4 .venv/bin/python my_script.py`. A 1D run with a small atom takes 1-5 s; with a large one (`fe607a`) about 20 s at 4 threads. When a script needs more than a minute, give that Bash call `timeout: 600000` and wait for it in the foreground; never send it to the background and poll its output.
2. Always launch with `verbose=0`. The console log goes to `data/<name>.log`. When `dispatch.x` fails, `tm.launch` raises `RuntimeError` with the error report of the log (an unknown namelist key comes with the closest valid key or the list of valid keys).
3. Print only the final numbers. Never print arrays, namelists or whole logs.
4. Never parse files under `data/.../save/` yourself; `tm.read` decodes them.
5. `tm.launch` must run in the OpenMULTI3D directory (it calls `./dispatch.x` there).
6. Leave tracked files unchanged: pass changes as `nml_settings`, never edit a namelist in `examples/` or `benchmarks/`. 1D models only: never run 3D atmospheres on this machine.
7. A run "with defaults" starts from the from-scratch dict below. The namelists in `benchmarks/` add an abundance table, `relative_abundance` and an opacity table, which change the result; copy them only when asked for those settings.
8. Finish every part of the task (numbers and files) before you reply; your reply ends the task. Check a figure by its file size, not by viewing it.

## Launch

```python
import numpy as np, astropy.units as u, tumult as tm

OUT = 'scratch/<topic>'      # every run name starts with your output prefix
# a tracked namelist with overrides; {group: {key: value}}
folder, ran = tm.launch(f'{OUT}/ba241', default='examples/sun_ba_1d.nml',
                        nml_settings={'atom_params': {'abundance': 2.41}},
                        threads=4, verbose=0)
# from scratch: NLTE departures only (no output spectrum)
folder, ran = tm.launch(f'{OUT}/o_sun', nml_settings={
    'io_params': {'datadir': 'data', 'do_trace': False},
    'atmos_params': {'atmos_file': 'input/atmos/sun.mod'},
    'atom_params': {'atom_file': 'input/atoms/atom.o41f'},
    'nlte_params': {'convlim': 1e-3},
    'rt_output_params': {'maxiter': 0}}, threads=4, verbose=0)
# from scratch: LTE synthesis with a line list and no model atom
folder, ran = tm.launch(f'{OUT}/lte_6162', nml_settings={
    'io_params': {'datadir': 'data', 'do_trace': False},
    'atmos_params': {'atmos_file': 'input/atmos/sun.mod'},
    'linelist_params': {'linelist_file': 'input/linelists/SIU_masterline_Bergemann2026.txt'},
    'spectrum_params': {'R_sampling': 5e5, 'aa_blue': 6160.0, 'aa_red': 6164.0}},
    threads=4, verbose=0)
run = tm.read(folder)
```

- `aa_blue`/`aa_red` are air wavelengths in Å.
- Equivalent widths and profiles need a resolved spectrum: set `spectrum_params` `R_sampling` (5e5) with `aa_blue`/`aa_red` around the lines, unless the namelist you start from already sets it. Without it a line-list run has no spectrum, and an NLTE run samples each line only on the atom's frequency grid (9 points per line in `fe607a`, 201 in `ba06`; the EW of Fe I 6151.6 comes out 79.7 mÅ instead of 50.4 mÅ).
- The run lands in `data/<name>/`. `ran` is `False` when a finished run of that name exists; it is then reused, not recomputed.
- `examples/sun_ba_1d.nml`: Sun (`input/atmos/sun.mod`), six-level Ba atom (`input/atoms/atom.ba06`, `abundance=2.11`), NLTE with `convlim=1d-2`. Every NLTE run also holds the LTE spectrum.
- 1D models: `input/atmos/*.mod` (MARCS; `sun.mod`, `metal_poor_tos.mod`, ...). Atoms: `input/atoms/atom.*` and `input/atoms/public/atom.*`.
- `atom_params: abundance` sets A(X) = log10(N_X/N_H) + 12 of the atom.
- Several runs: loop over `tm.launch` in the same script, each with its own name.

## Read

```python
line = run.get_line(4934.077)                 # model-atom transition nearest this air wavelength [Å]
w_nlte = line.crop().get_weq(target_unit=u.mAA).value
w_lte = line.crop(LTE=True).get_weq(target_unit=u.mAA).value
spec = run.crop(window=(4553.0, 4555.0))      # normalised flux in [lo, hi] Å; LTE=True for LTE
spec = run.crop(x0=4554.0, width=2.0)         # the same window by its centre and full width
x, y = spec.xdata.value, np.asarray(spec.ydata, float).ravel()   # Å, continuum = 1
w = np.trapezoid(1 - y, x) * 1e3              # equivalent width over the window [mÅ]

lvl = run.atoms[0].lvl                        # pandas: ev, g, ion (1 = neutral), label
ltau, b = run.get_depart(0)                   # level index from 0; ltau = log10 tau_500
ltau, b = np.asarray(ltau, float), np.asarray(b, float).ravel()
o = np.argsort(ltau)
logb = np.interp(-1.0, ltau[o], np.log10(b[o]))   # log10 b at log tau_500 = -1
grounds = [int(np.flatnonzero(lvl.ion.values == k)[0]) for k in np.unique(lvl.ion)]
```

`get_line` knows only transitions of the model atom; read lines of the line list with `run.crop(window=...)`. In `crop`, `dx` is the half-width and `width` the full width around `x0`. `ltau` runs from the bottom of the atmosphere upward, so sort it before `np.interp`.

## Abundance from an equivalent width

NLTE corrections of a few tenths of a dex are common, so span A0 - 0.5 ... A0 + 0.2 (for a large atom, first run A0 to see the sign of the correction, then 4-5 points in 0.1 dex steps on that side):

```python
A0, grid, w = 2.11, np.round(2.11 + np.arange(-0.5, 0.21, 0.05), 2), []
for a in grid:
    f, _ = tm.launch(f'{OUT}/ba_{a:.2f}', default='examples/sun_ba_1d.nml',
                     nml_settings={'atom_params': {'abundance': float(a)}}, threads=4, verbose=0)
    w.append(tm.read(f).get_line(5853.7).crop().get_weq(target_unit=u.mAA).value)
w = np.array(w)                               # W grows with A
assert w[0] < w_target < w[-1], 'widen the grid'
a_fit = np.interp(np.log10(w_target), np.log10(w), grid)
```

## Plot

Before writing figure code, read the `## Rules` section of `.claude/skills/pyplot-style/SKILL.md` (Read with a line limit, up to the next `##` heading) and follow every rule, including panel letters on multi-panel figures. Then write and run the script in one Bash call.

Get pyplot from tumult, which applies the project style, and give every figure the constrained layout:

```python
from tumult.plotting_scripts import get_plt
plt = get_plt()                               # pyplot in the project style
fig, ax = plt.subplots(layout='constrained')
```

A line profile, LTE against NLTE (`lte`, `nlte` from `line.crop(...)`):

```python
ax.plot(lte.xdata, lte.ydata, '--', color='tomato', label='1D LTE')
ax.plot(nlte.xdata, nlte.ydata, color='black', label='1D NLTE')
ax.set_xlabel('Wavelength [Å]')
ax.set_ylabel('Normalized Flux')
ax.legend(handlelength=1.5)
fig.savefig(f'{OUT}/profile.pdf')
```

Flux axes read exactly `Normalized Flux`, whatever the task calls the quantity. Departure plots: `ax.plot(ltau, np.log10(b))` per level, labelled with `lvl.label[i]`, x label `$\log_{10}\,\tau_{500}$`, y label `$\log_{10}\,b$`.

## Convergence

The NLTE iteration table is in the run log `data/<name>.log`, one row per iteration (`iter/maxiter`, `*` after an Ng step, largest relative change in %, its level, its cell):

```sh
grep -E '^ *[0-9]+\*? */ *[0-9]+ ' data/<name>.log | tail -3
```

`data/<name>/convergence_report_000.txt` lists the largest change of each level in the last iteration (levels numbered as in the table).

## On failure

- A `RuntimeError` from `tm.launch` quotes the error report of `dispatch.x`; fix the setting it names and rerun. Warnings stay in the log: `grep -m5 WARNING data/<name>.log`. `OpenMULTI3D terminated` in the log means the run finished.
- A traceback in your script: fix that line and rerun the same script; finished runs are reused, so a rerun costs only the reading. Edit the script with Edit or Write: `sed -i` on this machine (macOS) needs `sed -i ''`.
