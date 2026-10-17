---
name: tumult-plots
description: Read and plot OpenMULTI3D run output with tumult (spectra, line profiles, departure coefficients, rates, convergence, contribution functions, atmospheres, curves of growth, Grotrian diagrams). Use whenever plotting or inspecting the results of a dispatch.x run or a tumult workflow, in a script, notebook or manual figure.
---

# Plotting tumult output

The manual (`docs/manual/main.tex`) documents how to read and plot a run. Read the section
for the quantity before writing code, and use the maintained functions it names. Figure
style (colours, sizes, labels, panel letters) is in the `pyplot-style` skill.

| Quantity | Manual section (label) |
|---|---|
| reading a run, spectra, `crop` | Reading and Plotting Results (`sec:python`), Cropping options and Sets, modes and departures (`sec:pyref`) |
| plotting methods, overview table | Plotting (`sec:plotting`) |
| observed spectra, equivalent widths, fits | Observed spectra and measurements (`sec:measure`) |
| departure coefficients | Reading and Plotting Results; Saved populations (`sec:savedpops`) |
| NLTE convergence | NLTE convergence (`sec:convergence`) |
| contribution functions | Contribution Functions (`sec:cntrbf`) |
| centre-to-limb variation | Centre-to-limb variation |
| 3D atmosphere | Inspecting a snapshot |
| curves of growth, NLTE corrections | the worked example (`sec:example`) |
| workflow results (`SynthRequest`, `Departures`) | Python Workflows (`sec:pyworkflows`), `docs/FAQ.md` |

`docs/manual/figs/science_case.py` draws every figure of the manual from runs and is the
worked reference for combining these methods in one figure.

## Reading

- Read a run with `tm.read(folder)` and spectra with `run.crop(...)` or `line.crop(...)`.
  Never parse `save/flux_*.bin` or other `save/*.bin` files directly; arrays without a
  reader are assembled with `run.read_patch_save(label, zdim=...)`.
- Departure coefficients: `run.get_depart(level, atom)`; from a workflow, a solved
  `Departures` gives `.profile()` (`tumult/departures.py`).
- Rates: `tumult.atom_rates.read_rates(run_dir, atom_id)` reads the files written with
  `&output_params save_net_rates=T` (first column of the run only); `net_rates` and
  `transition_relevance` in the same module reduce them.
- Convergence: `run.conv_report` holds the `convergence_report_*.txt` tables (set when
  `run.atoms` is first read), and
  `Departures.convergence_levels()` the final correction per level. The per-iteration
  history is only in the console log (`rank_NNNNN.log`).

## Reuse first

Look for a plotting function made for the purpose before writing one, and extend it where
it falls short:

- `tumult.Run`: `plot_depart`, `plot_resolved`, `plot_cog_corr`, `plot_line_labels`,
  `plot_line_stem`, `plot_cntrbf`, `plot_cntrbf_map`
- `Spectrum.plot`, the atmosphere's `plot_slice`, `plot_ltau_mean`, `plot_ltau_cube`,
  `plot_geo_cube`, `plot_ltau_hist`, the transition's `plot_clv` and `plot_resolved`,
  the atom's `plot_grotrian`
- `tumult/plotting_scripts.py`: `comparison_scatter`, `scatter_density`, `quickplot`,
  `plot_line_annotations`, `plot_arrow`
- `docs/manual/figs/science_case.py` for the manual's figures

A new plot of a run quantity that will be needed again belongs in tumult as a method next
to these, with its reading logic in the matching `get_*` method, rather than in a script.
