# Installing OpenMULTI3D

OpenMULTI3D synthesises stellar spectra from 1D and 3D model atmospheres, in LTE and
NLTE. A working installation has three parts:

- **`dispatch.x`**, the native executable built from `SRC/` and `HOP/`. It solves the
  radiative transfer for one namelist and writes the result to a run folder under
  `data/`.
- **`tumult`**, the Python package that prepares and launches runs, reads their output
  and provides the analysis workflows.
- **TumultStudio**, the optional browser GUI (Bokeh) built on `tumult`.

This guide ends with a verified installation. [QUICKSTART.md](QUICKSTART.md) then
introduces the calculations and [docs/GUIWorkflow.md](docs/GUIWorkflow.md) the GUI.

## 1. Checkout layout and requirements

OpenMULTI3D is an experiment inside the DISPATCH repository and is built with
DISPATCH's shared build configuration. Clone the repository and work in the experiment
directory; every command in this guide runs there:

```bash
git clone https://github.com/Bergemann-Group-development/dispatch.git
cd dispatch/experiments/OpenMULTI3D
```

The local [Makefile](Makefile) includes `../../config/Makefile` from the parent
checkout, so a copy of this directory on its own does not build.

You need:

- macOS or Linux. The commands below were verified on macOS 15 (Apple silicon, Homebrew
  GCC 15); Linux x86-64 uses the same commands.
- A Fortran compiler: `gfortran`, `ifort` or `ifx`. On Apple silicon use gfortran: the
  macOS Intel compilers produce x86-64 code only, which runs slower under Rosetta, and an
  x86-64 HOP library cannot be loaded by an arm64 Python at all.
- `make` and `perl` (the build scans module dependencies with a Perl script).
- Python 3.11 or later (verified with 3.14), with `venv` and `pip`.

```bash
brew install gcc            # macOS: provides gfortran
sudo apt install gfortran   # Debian/Ubuntu
gfortran --version          # check it is on PATH
```

## 2. Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .          # headless use: the tumult package
pip install -e ".[gui]"   # adds TumultStudio (Bokeh) and notebook support
```

The editable install (`pyproject.toml`) makes `tumult` and `TumultStudio` importable from
any directory while the environment is active. They still read `libhop.so`, `dispatch.x`
and `input/` from this checkout, so do not move it.

## 3. Build the native code

```bash
make -j
```

The compiler is detected (`gfortran`, otherwise `ifort`); `make -j COMPILER=ifort`
forces one, and everything below is built with it. About two minutes on eight cores.
The build produces

- `dispatch.x` in this directory;
- `HOP/standalone/build/libhop.so`, the EOS and opacity library that `tumult.hop` loads
  through the committed header `HOP/standalone/libhop.h`;
- `HOP/standalone/build/hop_prune.x`, the line-list pruner the synthesis workflows run.

`make clean && make -j` rebuilds everything from scratch; `make -C HOP/standalone`
rebuilds the HOP library alone.

## 4. Input data

The repository tracks the inputs the verification below needs: one MARCS solar model,
four model atoms (Ba, Fe, Mg, O), abundance tables, continuum and molecular opacity data
and the atomic line list `input/linelists/master_line.txt`.

The larger inputs are not in git. They can be fetched from TumultStudio's Preferences
or from the command line, whichever suits; both use the same code and skip files that are
already complete, so an interrupted download resumes.

| Data | Size | TumultStudio page | Command line |
|---|---|---|---|
| Per-species master line lists | about 1 GB | Preferences → Line lists | `python -m tumult.synthesis.input_data fetch linelists` |
| Public model atoms (17 atoms) | 81 MB | Preferences → Model atoms | `python -m tumult.synthesis.input_data fetch atoms` |
| Metal-opacity (absmet) tables | four files of 182 MB | Preferences → Opacity tables | `python -m tumult.synthesis.input_data fetch absmet` |
| MARCS grid, standard composition | about 400 MB | Preferences → Atmosphere grids | `python -m tumult.synthesis.grid_catalog fetch marcs input/atmos/MARCS` |
| Stagger and M3DIS 3D grids | per model and snapshot | Preferences → Atmosphere grids | `python -m tumult.synthesis.grid_catalog fetch <grid> <dest> --snapshots N` |

The line lists, atoms and tables come from one MPCDF Keeper share
(<https://keeper.mpdl.mpg.de/d/3191475aebe340638422/>). Without `--dest`, `input_data`
puts each dataset where the code looks for it: `input/linelists/master_linelists/`,
`input/atoms/public/` and `input/absmet_files/`. `list` shows the files and sizes, `fetch`
takes `--file NAME` (repeatable) to download only some, and `--dry-run` prints the plan.
`python -m tumult.synthesis.grid_catalog list` names the grids; its module docstring
describes the model and snapshot selection. Downloaded files are git-ignored. The
verification below needs none of them.

In TumultStudio each page lists the share when opened and offers a download or delete
button per file, or for all files. The Line lists and Opacity tables pages download into
the folder their Synthesis preference (*Line list*, *Absmet*) names, and into the stock
folder above when that preference is a single file. The Model atoms page always uses
`input/atoms/public/`, where the atom picker of a calculation setup starts.

## 5. Verify the installation

Run the demo namelist:

```bash
OMP_STACKSIZE=256M ./dispatch.x demo.nml
```

`demo.nml` synthesises the four Ba II lines of a six-level toy model atom
(`input/atoms/atom.ba06`: 4554, 4934, 5853 and 6497 Å) in a MARCS solar model, in LTE
and NLTE, without a background line list. It takes well under a minute (about one second
on eight cores). `tm.launch` sets the same `OMP_STACKSIZE` for every run it starts.
The run succeeded when the last line of the output is

```
OpenMULTI3D terminated successfully.
```

and `data/demo/` holds the run. Then read the result from Python:

```python
import tumult as tm

run = tm.read("data/demo")
for lte in (True, False):
    core = run.crop(x0=4554.033, dx=1.0, LTE=lte).ydata.min()
    print("LTE " if lte else "NLTE", f"Ba II 4554 core flux {core:.3f}")
print(tm.HopEOS(feh=0.0).load_state([5777.0], rho=[2.5e-7]).ne)
```

Expected output (the first import compiles numerical kernels and takes a few seconds):

```
LTE  Ba II 4554 core flux 0.082
NLTE Ba II 4554 core flux 0.030
[2.02578179e+13]
```

The first two lines are the normalised flux in the line core; the NLTE core is deeper
than the LTE one for this atom. The last line is the electron density from the standalone
HOP library and shows that `libhop.so` loads. Your installation is verified when all
three lines print with these values.

## 6. Next steps

- [QUICKSTART.md](QUICKSTART.md) continues from this calculation: it plots the spectrum,
  compares LTE with NLTE, changes an input and opens the result in the GUI, with optional
  examples of the Python workflow interface and a 3D snapshot. [docs/FAQ.md](docs/FAQ.md)
  answers further workflow questions (departures, long jobs, abundances).
- TumultStudio, the GUI:

  ```bash
  python -m TumultStudio
  ```

  This serves the app behind a landing page at `http://localhost:5006/` and opens it in
  the browser; `--port N` changes the port and `--no-show` skips the browser.
  `bokeh serve --show TumultStudio` runs the app without the landing page. The GUI keeps
  its preferences, recent files and saved workspaces in `.tumultstudio.json` in the
  directory it was started from (`TUMULT_CONFIG` names another file); the file is
  git-ignored. [docs/GUIWorkflow.md](docs/GUIWorkflow.md) describes the tools.

## 7. Troubleshooting

**`make: *** No Fortran compiler found`**: neither `gfortran` nor `ifort` is on `PATH`.
Install one (section 1) or pass `COMPILER=<name>`; compiler profiles live in
`../../config/compiler/`.

**`make: *** No rule to make target '/Makefile.dep'`**: the parent DISPATCH checkout
is missing, so `../../config/Makefile` was not found. Build inside
`dispatch/experiments/OpenMULTI3D` of a full clone (section 1).

**`CERTIFICATE_VERIFY_FAILED: unable to get local issuer certificate`** from
when downloading input data on macOS: the python.org framework Python ships without a CA
trust store. Run the installer's certificate step once, for example
`/Applications/Python\ 3.12/Install\ Certificates.command`, or run with the virtual
environment active: `pip install -e .` installs the `certifi` bundle, which the downloads
use whenever it is importable.

**`ModuleNotFoundError: No module named 'numba'`** (or another installed package) when
running a script: a conda environment prefixed onto `PATH` shadows `.venv` with an
interpreter that carries none of the dependencies. Run `.venv/bin/python` explicitly, or
deactivate conda before activating `.venv`.

**`ModuleNotFoundError: No module named 'tumult'`**: the environment with the editable
install is not active. Activate `.venv` (section 2).

**`libhop.so not found at HOP/standalone/build/libhop.so`**: the standalone library was
not built, or was built into another directory. Run `make -C HOP/standalone`, or point
`HOP_LIB` at an existing library.

**`OSError: ... Address already in use`** when starting TumultStudio: another server
holds port 5006. Start with `--port N`.
