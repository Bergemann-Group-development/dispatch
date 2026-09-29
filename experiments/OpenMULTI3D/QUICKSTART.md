# OpenMULTI3D quickstart

This tutorial runs one small spectrum synthesis, reads its output in Python, plots it,
compares LTE with NLTE, changes an input and opens the result in TumultStudio. It needs
the verified installation from [README_INSTALL.md](README_INSTALL.md): the built
`dispatch.x`, the Python environment with `requirements-gui.txt` installed and the
tracked input data. No download is required until the optional 3D example.

## Where to work

Every command and code block runs in this directory (`dispatch/experiments/OpenMULTI3D`),
with the Python environment activated; all paths are relative to it. A script or notebook
saved here imports `tumult` directly; one kept elsewhere needs this directory on its
Python path, for example through `PYTHONPATH`. New files here and the outputs under
`data/` are ignored by Git, so following the tutorial leaves the tracked files unchanged.

Shell commands are shown with a `$` prompt; everything else is Python.

## 1. The example

`demo.nml` is the namelist the installation check runs. It names two inputs:

- **The atmosphere** `input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_..._s+0.00.mod`, the
  MARCS model of the Sun: a plane-parallel, hydrostatic 1D structure with
  Teff = 5777 K, log g = 4.4, [Fe/H] = 0 and a microturbulence of 1 km/s (`t01`), which
  the run reads from the file.
- **The model atom** `input/atoms/atom.ba06`, a six-level barium atom in MULTI format:
  the Ba I ground state, four Ba II levels (6s, 5d 2D3/2, 6p 2P1/2, 6p 2P3/2) and the
  Ba III ground state, coupled by four radiative lines, the Ba II resonance lines at
  4554.03 and 4934.08 Å and the subordinate lines at 5853.67 and 6496.90 Å, five
  photoionisation continua and tabulated electron and hydrogen collision rates.

The atom is a teaching example. Its 5d 2D5/2 level is missing, so the 6141 Å line is not
included, and the lines carry no hyperfine or isotopic splitting, which matters for the
strong 4554 Å line in the real Sun. Its purpose is to show the LTE and NLTE workflow;
abundances derived from it are not to be used in a publication. The larger atoms under
`input/atoms/` (Fe, Mg, O) are the ones the analysis workflows use.

`demo.nml` computes the statistical equilibrium of this atom in NLTE, with the LTE
solution as a by-product, and writes the emergent spectrum over the four lines only: with
no `spectrum_params` section, the wavelength grid is the atom's own line quadrature, 201
points per line. No background line list is included (the `linelist_params` line is
commented out), so the lines sit on a pure continuum.

Wavelengths in OpenMULTI3D output and in `tumult` are **air wavelengths in Ångström**,
the convention of most optical line lists and observations. The model atom stores level
energies, from which the code derives vacuum wavelengths; `tumult` converts every
wavelength it reads to standard air, whatever its range.

## 2. Run the calculation

```
$ OMP_STACKSIZE=256M ./dispatch.x demo.nml
```

`dispatch.x` takes one namelist and writes its run under `data/<namelist name>/`, here
`data/demo/`. The console shows the atmosphere and atom being read, then one
`Iteration n/99` block per iteration of the statistical equilibrium, each ending with the
largest relative population change (`global`), until `Converged` appears once that change
falls below `convlim` (0.01 in `demo.nml`). The run ends with

```
OpenMULTI3D terminated normally.
```

after about one second on eight cores. The installation check created the same run;
running it again overwrites it. `OMP_STACKSIZE` gives the OpenMP threads the stack the
solver needs; the Python launcher used in step 6 sets it for you.

Before execution, `Memory estimate rank ...` reports the expected peak of the memory
Dispatch registers on each MPI rank, in GiB (2³⁰ bytes). The equation-of-state setup,
the NLTE iteration and the output pass run one after the other, so the estimate is the
persistent allocation plus the largest of the three phases. The EOS phase grows with
the cells of the composition patches that run at the same time. In the RT phases each
busy thread sweeps one column (a frequency window and ray direction through all depth
slabs) and holds the workspace of one slab and the column's boundary buffer; with
`Jnu_to_file=F` the radiation fields of every slab are held as well. Raising `n_nu`
splits the frequencies over more columns and shrinks both per-thread terms.

At completion, `Memory comparison rank ...` gives the estimate, the registered peak and
their difference. Opacity caches, solver temporaries and allocations Dispatch does not
register are left out, and the process RSS is correspondingly higher.

## 3. Load the output

```python
import tumult

run = tumult.read("data/demo")
print(run.mode, run.lam.size, run.lam.min(), run.lam.max())
print(run.abundances["Ba"])
```

```
NLTE 804 4552.817856352975 6497.548103144807
2.11
```

`run` holds the run's namelist (`run.nml`), the composition it used (`run.abundances`),
the wavelength grid (`run.lam`), the emergent intensities and fluxes and the model atom
(`run.atoms[0]`). `run.mode` says whether the run solved NLTE or was an LTE run, which
carries no separate LTE spectrum. The abundance is on the
astronomical scale, A(Ba) = log10(N_Ba/N_H) + 12 = 2.11, taken from the atom file because
`demo.nml` sets `use_atom_abnd=T`; the standard composition
`input/abund/abund_magg22` has 2.17.

The 804 points are the four line windows, so the grid has gaps between them. Never read
the binary files under `data/demo/save/` yourself: `run` decodes their layout and unit
conventions.

## 4. Plot the normalised flux

```python
import matplotlib.pyplot as plt

spec = run.crop()          # NLTE flux divided by the continuum flux
spec.plot(color="black")
plt.show()
```

`run.crop()` returns a `Spectrum`, a wavelength and flux array pair carrying astropy
units. Two conventions apply:

- **Flux versus intensity.** The run computes the specific intensity I_λ along several
  outgoing directions μ = cos θ (`run.mus`, four values here) and integrates them over the
  stellar disc into the flux F_λ. `crop()` returns the flux by default; `imu=` selects one
  direction instead (optional example below). The physical units are
  erg s⁻¹ cm⁻² Å⁻¹, per steradian for an intensity; `run.get_spectrum()` returns the flux
  with those units.
- **Normalisation.** `norm=True`, the default, divides by the continuum flux computed in
  the same run at every wavelength, so the result is dimensionless with the continuum at
  1. Observations are normalised by fitting a continuum instead; the two are compared on
  that scale.

The plot shows the continuum at 1 with four absorption lines. Because the grid holds only
the four line windows, the segments between windows are straight lines drawn across the
gaps. To look at one line, crop around it:

```python
run.crop(x0=4554.033, dx=1.5).plot()
plt.show()
```

`x0`/`dx` crop to x0 ± dx; `xmin=`/`xmax=` or `window=(xmin, xmax)` are the alternatives.
`relative=True` puts the x axis in Δλ from `x0`, and `LTE=True` selects the LTE
spectrum.

## 5. LTE against NLTE for one line

The model atom's transitions are objects that know their own wavelength window:

```python
line = run.get_line(4554.033)          # the transition nearest that wavelength
nlte = line.crop()
lte = line.crop(LTE=True)

nlte.plot(label="NLTE")
lte.plot(linestyle="--", label="LTE")
plt.title(line.title)
plt.legend()
plt.show()

print(f"core flux  LTE {lte.ydata.min():.3f}  NLTE {nlte.ydata.min():.3f}")
```

```
core flux  LTE 0.082  NLTE 0.030
```

The NLTE core is deeper than the LTE core, the usual behaviour of a strong resonance
line: photons escape through the line core, the source function falls below the Planck
function in the upper photosphere where the core forms, and the line darkens. The wings
are nearly identical, since they form deeper in, where collisions keep the populations
close to LTE. The other three lines behave the same way; `run.atoms[0].line` lists them,
and a `Spectrum` measures its equivalent width with

```python
import astropy.units as u

for tr in run.atoms[0].line:
    w_lte = tr.crop(LTE=True).get_weq(target_unit=u.mAA)
    w_nlte = tr.crop().get_weq(target_unit=u.mAA)
    print(f"{tr.lam0:9.3f}  LTE {w_lte:.1f}  NLTE {w_nlte:.1f}")
```

```
 4554.033  LTE 138.7 mAngstrom  NLTE 142.6 mAngstrom
 4934.077  LTE 116.0 mAngstrom  NLTE 122.8 mAngstrom
 5853.675  LTE 59.1 mAngstrom  NLTE 62.8 mAngstrom
 6496.898  LTE 82.2 mAngstrom  NLTE 93.5 mAngstrom
```

The exact values depend on the compiler and platform in the last digit.

## 6. Change the barium abundance

`tumult.launch` writes a namelist from a default and a dictionary of overrides, runs
`dispatch.x` on it and returns the output folder. The tracked `demo.nml` stays as it is;
the new namelist is written under the name you give:

```python
folder, ran = tumult.launch("quickstart_ba_rich", default="demo.nml",
                            nml_settings={"atom_params": {"abundance": 2.41}},
                            verbose=1)
rich = tumult.read(folder)
print(rich.abundances["Ba"])
```

```
quickstart_ba_rich.nml launching
2.41
```

This writes `quickstart_ba_rich.nml` here, runs it with the console output in
`data/quickstart_ba_rich.log` and leaves the run in `data/quickstart_ba_rich/`. The
`atom_params` `abundance` entry overrides the atom file's abundance, so A(Ba) rises by
0.3 dex. `ran` is `False` when a finished run of that name already exists, in which case
nothing is recomputed. Any namelist entry can be overridden the same way, for example
`{"nlte_params": {"convlim": 1e-3}}` for a tighter convergence limit. Editing a copy of
`demo.nml` by hand and running the copy with `dispatch.x` is equivalent.

```python
line.crop().plot(label="A(Ba) = 2.11")
rich.get_line(4554.033).crop().plot(label="A(Ba) = 2.41")
plt.legend()
plt.show()
```

The richer line is deeper in the core (0.023 against 0.030) and wider in the wings; its
NLTE equivalent width grows from 143 to 179 mÅ, 0.10 dex for 0.30 dex in A(Ba). A line
this strong is saturated: its equivalent width grows more slowly than the abundance and
depends on the microturbulence and the damping as well, which is why an abundance
analysis prefers weaker lines or measures several of differing strength.

## 7. Save the result and open it in TumultStudio

A `Spectrum` can be written as a two-column text file that TumultStudio and other
tools read:

```python
import numpy as np

np.savetxt("quickstart_ba4554_nlte.txt", np.column_stack(nlte.values),
           header="air wavelength [AA]   normalised flux")
back = tumult.load_file("quickstart_ba4554_nlte.txt")
```

`Spectrum.save(path)` and `Spectrum.load(path)` round-trip the same object with its units
through JSON for use from Python. The run folder itself is the complete record: a later
session repeats step 3 on `data/demo` or `data/quickstart_ba_rich`.

To look at the runs in the GUI, start TumultStudio from this directory:

```
$ python -m TumultStudio
```

The landing page opens in the browser. Choose **Open ▾ → Spectrum…**, navigate the file
picker to this checkout's `data/demo` folder and press **Open run**: the run's normalised
NLTE flux appears on the canvas and in the object list on the left. Open
`data/quickstart_ba_rich` the same way to overlay the abundance-changed run, and open the
saved text file through the same picker with **Open** to add the cropped line. The canvas
tools zoom, measure and compare the spectra; [docs/GUIWorkflow.md](docs/GUIWorkflow.md)
describes them and the synthesis and analysis tools that run the same calculations from
the GUI.

**This is the end of the quickstart.** To remove what it created:

```
$ rm -rf quickstart_* data/quickstart_* data/demo
```

Where to go next:

- The optional examples below cover viewing angles and broadening, the Python workflow
  interface and a 3D snapshot.
- [docs/GUIWorkflow.md](docs/GUIWorkflow.md) for the GUI tools.
- [docs/FAQ.md](docs/FAQ.md) for scripting departures, long calculations, abundance and
  stellar-parameter analyses and the run cache.
- [docs/GLOSSARY.md](docs/GLOSSARY.md) for the terms these documents share.

---

## Optional: intensities at different viewing angles and broadening

The flux in step 4 integrates intensities along the directions listed in `run.mus`. Each
can be plotted alone:

```python
line = run.get_line(4554.033)
for imu, mu in enumerate(run.mus):
    line.crop(imu=imu, dx=0.8).plot(label=f"μ = {mu:.2f}")
plt.title(line.title)
plt.legend()
plt.show()
```

The normalised profiles of this line differ little between the directions: the core is
0.026 at disc centre and 0.047 at μ = 0.21, where the inclined ray sees the continuum
form higher up, closer to the layers that form the core. `norm=False` shows the
intensities themselves, with the continuum darkening towards the limb. The four
directions are the run's Lobatto quadrature (`rt_output_params`); a 1D model sees them
as four inclined rays through the same stratification, and the centre-to-limb variation
becomes a sharper test of a model in 3D, where the rays cross different granules.

A synthetic spectrum is sharper than any observed one. `crop()` convolves the extracted
spectrum before returning it: `vdop=` with a Gaussian of that Doppler velocity in km/s,
`res=` with the Gaussian of a spectrograph's resolving power, `vrot=` with a rotation
profile, and `vshift=` applies a radial-velocity shift:

```python
for vdop in [1.0, 3.0, 5.0]:
    line.crop(dx=0.5, vdop=vdop).plot(label=f"v = {vdop:.1f} km/s")
plt.legend()
plt.show()

line.crop(res=20_000).plot(label="R = 20 000")
```

The equivalent width is preserved by a convolution; the core depth is not, which is why
the profile comparison in step 5 is only meaningful at the same broadening.

## Optional: a synthesis through the Python workflow interface

Steps 2 and 6 drove `dispatch.x` through namelists. The analysis workflows (`tumult.synthesis`,
`tumult.departures`, `tumult.abundance`, `tumult.sparams`) and the GUI describe a
calculation as a `SynthRequest` instead: a configuration, an atmosphere, stellar
parameters, wavelength windows and, for NLTE, the model atoms to solve. The layer writes
the namelists, prunes the line list to the windows, caches every native run by its inputs
and returns `Spectrum` objects.

This example synthesises the same Ba II 4554 Å region on the same MARCS model, now
with the blends of the tracked line list `input/linelists/master_line.txt`:

```python
from tumult import departures, synthesis
from tumult.synthesis import Atmosphere, Star, SynthConfig, SynthRequest

config = SynthConfig(root=".", cache="quickstart_synth",
                     linelist="input/linelists/master_line.txt",
                     abundances="input/abund/abund_magg22",
                     resolving_power=200_000)
marcs = Atmosphere("input/atmos/p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_"
                   "n+0.00_o+0.00_r+0.00_s+0.00.mod")
request = SynthRequest(config=config, atmosphere=marcs,
                       star=Star(teff=5777, logg=4.44, feh=0.0, vmic=1.0),
                       windows=(4550, 4558))
request.check()                       # names any missing input before running

lte = synthesis.synthesize(request)
nlte_request = request.update(
    departures=departures.request(request, atom="input/atoms/atom.ba06"))
nlte = synthesis.synthesize(nlte_request)

lte.spectrum.plot(label="LTE")
nlte.spectrum.plot(label="NLTE")
plt.legend()
plt.show()
print(lte.runs, nlte.provenance["departures"])
```

The stellar parameters are always given explicitly, since the workflow never reads them
from a model file; `abund_magg22` is the standard composition, so Ba is at A(Ba) = 2.17
here rather than the atom file's 2.11. `lte` and `nlte` are `SynthesisResult` objects:
`.spectrum` is the normalised flux over the window, `.runs` the native run folders under
`data/quickstart_synth/` (readable with `tumult.read`, as in step 3) and `.provenance`
records the line-list pruning and whether departures were solved or reused. The first
call takes a few seconds; repeating either call returns the cached result at once. The
NLTE request solves the barium departures once and applies them in the formal solution
that produces the spectrum; the LTE request skips the solve.

An abundance change is an `abundances=` entry of the request, in [X/Fe] by default:
`request.update(abundances={"Ba": 0.3})`. The same request runs unchanged on a MARCS grid
interpolated at the star (`AtmosphereGrid("input/atmos/MARCS")`, after
`fetch_input_data.py`), on a 3D snapshot (`Snapshot`), and as the input of the
abundance and stellar-parameter workflows; [docs/FAQ.md](docs/FAQ.md) has the
recipes. `rm -rf quickstart_synth data/quickstart_synth` removes this example's cache.

## Optional: a 3D snapshot

This example downloads 82 MB and runs for several minutes on eight cores.

`demo3d.nml` runs the same barium atom in one snapshot of a solar 3D radiation-hydrodynamic
simulation, in M3DIS format. Download and unpack it into `input/atmos/`, where downloaded
files are ignored by Git:

```python
import zipfile
from tumult.synthesis.remote import fetch_file

url = ("https://gemini-web.mpia.de/atmos3D/SunLTE_E_t57.77g44.40m0.000_v1.0/"
       "SunLTE_t5777_g4.440_z0.000_sn779_m3dis.zip")
archive = fetch_file(url, "input/atmos/SunLTE_t5777_g4.440_z0.000_sn779_m3dis.zip")
with zipfile.ZipFile(archive) as zf:
    zf.extractall(archive.parent)          # the .atmos cube (98 MB) and its .mesh
archive.unlink()
```

`fetch_file` is the downloader behind the grid and line-list catalogues: it times out
on a stalled connection instead of hanging and resumes a partial download when the
block is run again.

The snapshot is a 160 × 160 × 159 cube; `demo3d.nml` samples it onto 30 × 30 columns
(`atmos_params` `nx`, `ny`) of 96 cells (the `nz` default), coarse enough for a
demonstration. Run it:

```
$ OMP_STACKSIZE=256M ./dispatch.x demo3d.nml
```

The iteration blocks look like the 1D ones, and each takes far longer: every iteration
solves the transfer along the angle quadrature through 900 columns instead of one, and
the solution converged after 13 iterations in six and a half minutes on an eight-core
laptop. The output is `data/demo3d/`. `tumult.read` opens it like a 1D run; the saved
intensities and the flux are means over the horizontal plane, so the profile is the
disc-integrated line of the whole snapshot:

```python
run1d = tumult.read("data/demo")
run3d = tumult.read("data/demo3d")

run1d.get_line(4554.033).crop().plot(label="1D MARCS, NLTE")
run3d.get_line(4554.033).crop().plot(label="3D snapshot, NLTE")
plt.legend()
plt.show()
```

The 3D core is much shallower (0.115 against 0.030) and the profile slightly asymmetric
and blueshifted, with a marginally smaller equivalent width (133 against 143 mÅ). The
convective velocity field replaces the 1 km/s microturbulence of the MARCS model, and
the up- and downflows Doppler-shift the contributions of the individual columns, so the
average of many shifted profiles has a broader, shallower core; the granulation's
temperature inhomogeneities change the line formation as well. A single snapshot at
this sampling is not a converged 3D result: a 3D abundance analysis averages over
several snapshots and checks the sampled mesh. The `Snapshot` atmosphere of the workflow
interface takes the mesh as its `resolution`, a snapshot grid averages a node's spectrum
over its selected snapshots (`tumult.synthesis.snapshot_grid`), and
`tumult.synthesis.grid_catalog` fetches the published 3D grids per model and snapshot.
