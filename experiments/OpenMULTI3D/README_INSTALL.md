# OpenMULTI3D — Getting Started

You've just cloned the repo. This experiment has two parts, both needed for a
full setup:

1. **Fortran radiative-transfer code** (`SRC/`, `HOP/`) — compiles into the
   `dispatch.x` executable that runs the synthesis.
2. **Python stack** — the `tumult/` synthesis library and the `TumultStudio2/`
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

This installs the headless `tumult/` stack. For TumultStudio2 and notebook
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

### Replaying a resolved synthesis

The higher-level API uses the same request for blocking and durable asynchronous work. A
request owns a copy of its masks and settings, so later GUI or script edits cannot alter a
submitted calculation:

```python
from tumult import GridTarget, MaskWindow, SynthRequest, open_workflow, submit

request = SynthRequest(
    prefs={"dispatch_cwd": ".", "synth_output_dir": "synth", "R": 300_000,
           "linelist": "input/linelist/vald_3000-10000.txt",
           "linelist_crop_aa": 2.0},
    run_vals={"mode": "NLTE", "atoms": ["input/atoms/MgI.atom"]},
    selection=GridTarget(
        {"folder": "input/atmos/MARCS",
         "pattern": "{Teff:d}_g{logg:+.1f}_m{:3}_t{vmic:02d}_st_z{FeH:+.2f}",
         "format": "Marcs"},
        (5777, 4.44, 0.0, 1.0)),
    masks=[MaskWindow(fit=(5167.0, 5168.0), synth=(5166.5, 5168.5))],
)

request.resolve().to_json("resolved-synthesis.json")  # inspect/replay exact inputs
job = submit(request)
print(job.id, job.status, job.progress)
result = open_workflow(job.folder).result()             # also works after reconnecting
result.spectrum.plot()
```

A grid of 3D models registers like a MARCS grid and interpolates between its nodes the same
way. Grid patterns match each file's path relative to the grid folder, listed recursively,
so a directory-per-node layout binds its parameters from the directory name;
`tumult.synthesis.STAGGER_PATTERN` covers the post-processed Stagger naming (`t50g47m00`,
`t5777g44m30`, `_new` models, `atm3d.*` snapshots). A pattern reads like the Python format
that wrote the name, `atm3d.t{Teff/100:02d}g{logg*10:02d}{FeH*10:+03d}{}_{snap:05d}`, with
head and tail implicitly wild and `p`/`m` accepted as signs; `compile_pattern` documents
the details (`|` alternatives, `{snap}`). A snapshot grid's record also
names the snapshots representing each node and, optionally, the mesh the cubes are
sampled onto:

```python
from tumult.synthesis import STAGGER_PATTERN

stagger = {"folder": "~/Atmospheres/Stagger/dwarfs", "format": "Stagger2",
           "pattern": STAGGER_PATTERN,
           "snapshots": "3",             # first three per node; "all"; or numbers "8,31"
           "resolution": [80, 80, 240]}  # blank: the nx/ny/nz preferences
marcs = {"folder": "input/atmos/MARCS", "format": "Marcs"}
selection = GridTarget(stagger, (5250, 4.5, 0.0, 1.0), preparation_grid=marcs)
```

Every selected snapshot becomes one 3D run; a node's spectrum is the mean over its snapshots
and the nodes combine with the usual interpolation weights. When several subfolders map to
one node, the lexically last wins (`t45g40m00_new` over `t45g40m00`). HOP line-list
preparation still needs 1D structures, hence the MARCS `preparation_grid`. Snapshot nodes
carry no vmic.

For a single 3D calculation use `Snapshot3D(path, stellar, resolution)`. Call
`tumult.synthesis.solve_departures(...)` without masks to create a `DepartureField`, save or
load it with `DepartureField.save/load`, and pass it as `departure_field=` when constructing
another request. The result records both the departure solve and formal-solution settings.

---

## 5. Launch TumultStudio2 (optional GUI)

From the repository root (the directory containing this file):

```bash
bokeh serve --show TumultStudio2
```

Run Bokeh from here so the repository package is importable.

The GUI keeps its preferences, recent files and saved sessions in `.tumultstudio.json`
in the directory Bokeh was started from (`TUMULT_CONFIG` names another file). The file
is created on first use; keep it out of version control.

### GUI tests

`python -m pytest tests/test_ui_execution.py tests/test_ui_events.py tests/test_store.py
tests/test_session.py` exercises the document, the job table, the config store and the
save slots headlessly, with stub workers and a temporary config. The browser half of
`tests/test_ui_events.py` additionally needs `node` and the `playwright` npm package with
its Chromium: it looks in the `npx` cache (`~/.npm/_npx/*/node_modules`) or under
`TUMULT_NODE_PATH`, and skips when neither has it. `npx playwright install chromium`
populates the cache; no Python Playwright is required.
