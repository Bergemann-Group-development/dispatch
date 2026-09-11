# TumultStudio2 workflow: ranges, line masks, synthesis and fits

Target design for the canvas tools and the synthesis entry points. It replaces the
Synthesize sidebar panel (`TumultStudio2/synthpanel.py`) and the separate EW toolbar tool
with one range gesture and one synthesis modal. Terms follow `GLOSSARY.md`.

## Left sidebar: one focused object

Spectra, line masks and results form one list with a single focus. Clicking an object
opens its view: a spectrum shows the inspector, a line mask the LineMasks view, a result
its table. No action relies on a second, implicit selection. Whenever an action needs a
spectrum in addition to a mask set or range (synthesis window, measurement reference, fit
target), it asks for it inside its modal, as the abundance and ion-balance modals do.

## Toolbar

`Annot | Filter [Filtered] Normalise | Mask | Corr`. The Synth and EW buttons are removed;
the Mask tool takes over both range gestures.

## Mask tool

1. Clicking **Mask** makes width-only box-select the active drag and the right sidebar
   shows "Drag a wavelength range". Esc or clicking **Mask** again leaves the tool and
   clears the band.
2. After a drag the band stays shaded and the sidebar shows the range as editable
   λ min / λ max fields with these actions:
   - **New line mask**: creates a blank line-mask set holding this window (fit window =
     range, synth window = range ± padding), adds it as a new object, focuses it and opens
     the LineMasks view, where the element is set. It never adds to an existing set.
   - **Synthesise…**: opens the synthesis modal with the window prefilled.
   - **Fit…**: opens the fit modal with the range as the fit window.
   - **Equivalent width**: the EW panel (profile selector, W, reduced W, hatched area)
     appears inline below the actions. It measures the focused spectrum and is disabled
     with a hint when the focused object is not a spectrum.
3. The band persists until the tool is left, so several actions can use one range.

## LineMasks view

Opened by focusing a line-mask object. The sidebar holds name and modified badge,
element, synth padding, save-as / revert, and the analysis actions **Synthesise
windows…**, **Fit windows…**, **Measure abundances** and **Fit ionisation balance**.
The centre shows the mask table below the spectrum; existing windows are resized by their
edge handles only.

**+ New mask** in the table arms the Mask tool's gesture in *add* mode: the sidebar shows
the drawn range with a single action, **Add window to <set>**. Leaving the tool returns to
the LineMasks view. Creating and extending a set therefore use the same drag.

## Synthesis modal

One component, opened from the Mask tool sidebar (free range), the LineMasks view (the
set's synth windows) and the empty landing page (blank λ fields).

- Header: what will be computed, e.g. "5700.00–5720.00 Å" or "7 windows of Fe-lines.txt,
  5480–6720 Å".
- Body: Atmosphere (grid, interpolate or pick, file), Stellar parameters, Abundances,
  Radiative transfer (model atoms in NLTE). Values persist in the `synth_run` store
  section so the modal reopens with the last configuration.
- Footer: **Synthesise**, Cancel. The result is a synth object on the canvas with the
  Info-tab provenance; a run from a mask set records the set name.

Departure-coefficient solves are not started from this modal; their entry point is decided
separately. Grid registration, paths and the run cache stay in Preferences → Synthesis.

## Fit modal

The same component in fit mode.

- Header: the fit window (range or windows) and a **Target spectrum** selector listing the
  loaded observations.
- Body: the four synthesis sections as the initial guess, plus a **Fit** section with the
  switches for Teff / log g / [Fe/H] / vmic (hidden for a fixed atmosphere file), the
  per-element abundance rows and vshift / vmac / vrot.
- Footer: **Fit**, Cancel. On completion the best-fit values are written back into the
  `synth_run` section, so the next synthesis starts from them; the fitted spectrum is added
  as a synth object and the broadening is applied as the live convolution in the Edit tab.

## Code consequences

- `SynthPanel` is replaced by a `SynthModal` structured like the abundance modal; `FitTab`
  becomes its Fit section; the departure-solve action is dropped from the GUI for now.
- `RangeSelect` is owned by the Mask tool alone. `EquivalentWidth` becomes a sub-panel of
  the Mask tool instead of a toolbar tool.
- `MaskTool._start_placing` is replaced by activating the Mask tool in add mode.
- The right sidebar gains a `mask-active` mode class; `synth-active` and `ew-active` go.
- The ion-balance and abundance modals are unchanged.
