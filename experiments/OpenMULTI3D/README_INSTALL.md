# OpenMULTI3D — Getting Started

You've just cloned the repo. This experiment has two parts, both needed for a
full setup:

1. **Fortran radiative-transfer code** (`SRC/`, `HOP/`) — compiles into the
   `dispatch.x` executable that runs the synthesis.
2. **Python stack** — the `tumult/` synthesis library and the `TumultStudio/`
   Bokeh GUI built on top of it.

Follow the steps below in order. Once set up, [docs/FAQ.md](docs/FAQ.md) answers
workflow questions (masked synthesis, 3D snapshots, departures, long jobs, abundances).

---

## 0. Prerequisites

- **Python 3.11+** (3.11 is the reference version; 3.12+ also works)
- **A Fortran compiler**: `gfortran` (recommended) or Intel `ifort`
- **git** (to have gotten this far) and **make**

Install a Fortran compiler if you don't have one:

```bash
# macOS (Homebrew)
brew install gcc          # provides gfortran

# Linux (Debian/Ubuntu)
sudo apt install gfortran
```

Check it's on your `PATH`:

```bash
gfortran --version
```

---

## 1. Set up the Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

This installs the headless `tumult/` stack. For TumultStudio and notebook
support, install the GUI extras instead:

```bash
pip install -r requirements-gui.txt
```

---

## 2. Compile the Fortran code

The build uses the top-level DISPATCH `config/Makefile` (two levels up)
together with the local [`Makefile`](Makefile). The compiler is auto-detected
(`gfortran` if present, otherwise `ifort`) — no manual configuration needed
on a standard macOS or Linux machine:

```bash
make -j
```

This produces the executable `dispatch.x` in this directory. To clean and
rebuild from scratch:

```bash
make clean
make -j
```

> If you have multiple Fortran compilers installed and want to force a
> specific one, pass it explicitly, e.g. `make -j COMPILER=ifort`. Compiler
> profiles live under `../../config/compiler/<name>/`.

---

## 3. Get the input data

Model atoms, line lists, MARCS atmospheres, and opacity/`absmet` files
(about 1.9 GB) are too large for git and are hosted on MPCDF Keeper:

<https://keeper.mpdl.mpg.de/d/3191475aebe340638422/>

The share mirrors the layout of [`input/`](input/). Download and merge it
into `input/` with the provided script (stdlib only, no extra packages
needed — works even before you set up the venv):

```bash
python fetch_input_data.py            # download + unpack + merge into ./input
python fetch_input_data.py --dry-run  # preview what would be downloaded
```

The script is re-runnable: files already present with the right size are
skipped, and archives are only unpacked once. The downloaded data stays
git-ignored, so it will not bloat the repository.

> **macOS TLS error** (`CERTIFICATE_VERIFY_FAILED: unable to get local issuer
> certificate`): the python.org framework Python ships without a CA trust store.
> Fix it once by running the installer's cert step, e.g.
> `/Applications/Python\ 3.12/Install\ Certificates.command`. The script also
> uses the `certifi` bundle automatically when it is installed, and offers a
> `--insecure` last resort that skips certificate verification.

---

## 4. Run a test synthesis

From this directory, with the input data in place:

```bash
OMP_STACKSIZE=256M ./dispatch.x demo.nml
```

`demo.nml` (or any of the other `*.nml` namelist files here) configures the
run. Output is written to `data/<input-path-without-suffix>/`; spectra land
under `save/` as `out_lam.bin` (the wavelength grid) plus `ie_0NNNNN.bin`
chunks.

See `demo.ipynb` for an example that runs `dispatch.x` and then loads and
plots the resulting spectrum with `tumult.read(...)`.

### Launching a run from Python

Instead of calling `dispatch.x` directly, `tumult.launch` can build a namelist
from a default `.nml` file plus overrides, run `dispatch.x` for you, and skips
the run entirely if a completed output already exists at that path:

```python
import tumult

output_folder, ran = tumult.launch(
    'demo_run',                                  # run name -> data/demo_run/, ./demo_run.nml
    default='demo.nml',                           # base namelist to start from
    nml_settings={'atom_params': {'convlim': 1e-3}},  # overrides on top of the default
)

run = tumult.read(output_folder)
run.get_spectrum(imu=0, LTE=True).plot()
```

`launch()` returns `(output_folder, ran)`, where `ran` is `False` if it served
a cached completed run instead of re-running `dispatch.x`. See
`tumult/launch.py` for the full set of options (`threads`, `run_async`,
`slurm_script`, etc.).

### Scripting a synthesis

A `SynthRequest` states a calculation in physical terms and is immutable, so later GUI or
script edits cannot alter a submitted calculation. The same request runs blocking or as a
durable workflow in the per-user daemon:

```python
from tumult import departures, synthesis
from tumult.synthesis import (AtmosphereGrid, Star, SynthConfig, SynthRequest,
                              open_workflow, submit_operation)

config = SynthConfig(root=".", cache="synth", resolving_power=300_000,
                     linelist="input/linelists/sun_3000-9000.txt",
                     abundances="input/abund/abund_magg22", linelist_skirt=2.0)
request = SynthRequest(config=config, atmosphere=AtmosphereGrid("input/atmos/MARCS"),
                       star=Star(5777, 4.44, 0.0, 1.0), windows=(5166.5, 5168.5))
request = request.update(departures=departures.request(request, atom="input/atoms/atom.mg86g"))

request.to_json("synthesis-request.json")          # SynthRequest.from_json reopens it
result = synthesis.synthesize(request)             # blocking
job = submit_operation("synthesis.synthesize", request)
result = open_workflow(job.folder).result()        # also works after reconnecting
result.spectrum.plot()
```

### Published grids

`tumult.synthesis.grid_catalog` fetches the published grids OpenMULTI3D reads and returns
the record that registers each: the standard-composition MARCS grid
(`marcs`, marcs.astro.uu.se), the Stagger 2.0 dwarfs and giants (`stagger-dwarfs`,
`stagger-giants`, Zenodo) and the M3DIS PLATO grid (`m3dis-plato`, gemini-web.mpia.de).
The 3D grids are fetched in parts, out of their archives, so a stellar parameter box and a
snapshot count per model keep a download small (one snapshot is 33 MB for Stagger and
100 to 200 MB for M3DIS; a whole grid at one snapshot per model is a few GB for Stagger and
about 30 GB for M3DIS):

```bash
bin/tpy -m tumult.synthesis.grid_catalog list
bin/tpy -m tumult.synthesis.grid_catalog models m3dis-plato
bin/tpy -m tumult.synthesis.grid_catalog plan  m3dis-plato input/atmos/M3DIS_PLATO \
    --snapshots 2 --teff 5000:6500 --logg 4:5
bin/tpy -m tumult.synthesis.grid_catalog fetch m3dis-plato input/atmos/M3DIS_PLATO \
    --model PLATO_t57.77g44.40m0.000_v1.0 --snapshots 2
```

`models` lists a grid; `plan` shows which models and snapshots a fetch would take without
transferring anything; `fetch` is resumable and records what it fetched in the folder's
`grid_manifest.json`. Running it again with a larger count only adds snapshots: a count
selects the leading entries of a fixed scrambled Halton ordering of each model's snapshots
(`tumult.synthesis.snapshot_grid.halton_order`), nested and spread over the time series,
and the run-time policy `"snapshots": "k"` selects the same entries. Grids fetched under
`input/atmos/` register themselves in TumultStudio. Its Preferences → Atmosphere grids
page sets one storage folder and one snapshot count for every grid, lists the supported
grids with what is on disk (download all, models, delete), and per grid the models with a
download or delete button each; custom grid folders register there too.

### Grid records

A grid of 3D models registers like a MARCS grid and interpolates between its nodes the same
way. Grid patterns match each file's path relative to the grid folder, listed recursively,
so a directory-per-node layout binds its parameters from the directory name;
`tumult.synthesis.STAGGER_PATTERN` covers the post-processed Stagger naming (`t50g47m00`,
`t5777g44m30`, `_new` models, `atm3d.*` snapshots) and `M3DIS_PATTERN` the M3DIS server
layout (`PLATO_t57.77g44.40m0.000_v1.0/m3dis_318.atmos`). A pattern reads like the Python
format that wrote the name, `atm3d.t{Teff/100:02d}g{logg*10:02d}{FeH*10:+03d}{}_{snap:05d}`,
with head and tail implicitly wild and `p`/`m` accepted as signs; `compile_pattern`
documents the details (`|` alternatives, `{snap}`). A snapshot grid's record also
names the snapshots representing each node and, optionally, the mesh the cubes are
sampled onto:

```python
from tumult.synthesis import STAGGER_PATTERN, AtmosphereGrid

stagger = AtmosphereGrid("~/Atmospheres/Stagger/dwarfs", format="Stagger2",
                         pattern=STAGGER_PATTERN,
                         snapshots="3",             # three per node in Halton order; "all"; or "8,31"
                         resolution=(80, 80, 240),  # None: the nx/ny/nz preferences
                         preparation=AtmosphereGrid("input/atmos/MARCS"))
```

Every selected snapshot becomes one 3D run; a node's spectrum is the mean over its snapshots
and the nodes combine with the usual interpolation weights. When several subfolders map to
one node, the lexically last wins (`t45g40m00_new` over `t45g40m00`). HOP line-list
preparation still needs 1D structures, hence the MARCS `preparation` grid. Snapshot nodes
carry no vmic.

For a single 3D calculation use `Snapshot(path, format, resolution)` in a `SynthRequest`.
`tumult.departures.solve(request, atom=...)` creates `Departures`, which reopen with
`departures.load` and are attached with `request.update(departures=...)`. The result
records both the departure solve and formal-solution settings. See
[docs/FAQ.md](docs/FAQ.md), "Scripting the workflows".

---

## 5. Launch TumultStudio (optional GUI)

From the repository root (the directory containing this file):

```bash
bokeh serve --show TumultStudio
```

Run Bokeh from here so the repository package is importable.

The GUI keeps its preferences, recent files and saved sessions in `.tumultstudio.json`
in the directory Bokeh was started from (`TUMULT_CONFIG` names another file). The file
is created on first use; keep it out of version control.

## 6. Tests

```bash
pip install -r requirements-test.txt
python -m pytest                    # headless gate: science, I/O and workflows
python -m pytest tests/studio       # GUI components; needs requirements-gui.txt
python -m pytest tests/browser      # browser journeys through TumultStudio
```

Run these from this directory. `pytest.ini` selects the headless directories by default;
every other suite is given by path. The browser suite needs `node` and the `playwright`
npm package with its Chromium: `npx playwright install chromium` populates the `npx`
cache (`~/.npm/_npx/*/node_modules`), and `TUMULT_NODE_PATH` names another location. No
Python Playwright package is required.

[tests/README.md](tests/README.md) describes the suites, their prerequisites and the
policy for admitting a test or a reference dataset.
