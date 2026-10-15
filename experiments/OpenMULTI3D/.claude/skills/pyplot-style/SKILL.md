---
name: pyplot-style
description: House style for matplotlib figures in OpenMULTI3D (notebooks, analysis scripts, docs/manual figures, papers and plotting code in tumult). Use whenever making or editing a plot with matplotlib.
---

# Matplotlib figures

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

## Rules

1. **Style:** `plt.style.use("tumult.Bergemann2020")` before the first figure.
   `tumult.plotting_scripts.get_plt()` resets matplotlib to its defaults and applies the
   same style on its first call. Call it before any `rcParams` change when the figure also
   uses tumult's plotting methods, otherwise their first call undoes the change.
2. **Colours:** the default colour cycle, in its order. Lines that may overlap differ in
   their dash style (`"-"`, `"--"`, `":"`, `"-."`) as well.
3. **Reference values** (observations, `axhline`/`axvline` guides, zero lines) are
   `color="lightgrey"`.
4. **Sizes:** the figure size sets the print size, and for LaTeX documents also the
   `\pltscale` the figure is included with. Change the font size where that is not enough.
   Leave tick and axis line sizes alone. Keep the standard `linewidth`; change it only as a
   last resort.
5. **Markers:** pick the size with the marker (`","`, `"."`, `"o"`) instead of
   `markersize` or `s`. A scatter of 3D data with one point per spatial grid point uses
   `marker=","` and `rasterized=True`.
6. **Logarithms:** an axis label of a logarithmic quantity says which one,
   `$\log_{10}$` or `$\ln$`, never a bare `$\log$`.
7. **Legends:** `handlelength=1.5`.
8. **Suptitles:** `fontsize="medium"`, `fontweight="bold"`.
9. **Output:** save as `.pdf`.

```python
import matplotlib.pyplot as plt

plt.style.use("tumult.Bergemann2020")
fig, ax = plt.subplots()
ax.plot(obs.x, obs.y, color="lightgrey", label="KPNO")
ax.plot(lte.x, lte.y, label="LTE")
ax.plot(nlte.x, nlte.y, "--", label="NLTE")
ax.axhline(1.0, color="lightgrey")
ax.set_ylabel(r"$\log_{10}\,(W/\lambda)$")
ax.legend(handlelength=1.5)
fig.suptitle("Mg I 5711", fontsize="medium", fontweight="bold")
fig.savefig(out / "mg5711.pdf")
```

## Figures for LaTeX documents

A LaTeX document includes every figure with `\includegraphics[scale=\pltscale]{...}`,
never `width=`, and has one `\pltscale`. The figure size is the printed size divided by
`\pltscale`, at full or half text width, with 4:3 panels. A row of three panels shares
its y axis. Save with `savefig.bbox` set to `"standard"` (the style sets `"tight"`), so the
PDF is exactly the figure size.

| Document | `\pltscale` | Plot text | Full width | Half width |
|---|---|---|---|---|
| `docs/manual` (16 cm text) | 1 | 8 pt | (6.30, 2.36) in | (3.15, 2.36) in |
| A&A paper (184 mm text) | 0.7 | style default | (10, 3.75) in per row | (5, 3.75) in |

Flux axes read "Normalized Flux". Residual or difference panels only on request.
