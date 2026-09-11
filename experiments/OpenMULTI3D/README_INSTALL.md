# OpenMULTI3D — Getting Started

You've just cloned the repo. This experiment has two parts, both needed for a
full setup:

1. **Fortran radiative-transfer code** (`SRC/`, `HOP/`) — compiles into the
   `dispatch.x` executable that runs the synthesis.
2. **Python stack** — the `tumult/` synthesis library and the `TumultStudio2/`
   Bokeh GUI built on top of it.

Follow the steps below in order.

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

For a 3D calculation use `Snapshot3D(path, stellar, resolution)`. Call
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
