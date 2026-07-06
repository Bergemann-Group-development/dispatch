# OpenMULTI3D — Installation & Build

This experiment has two parts:

1. **Fortran radiative-transfer code** (`SRC/`, `HOP/`) that builds into `dispatch.x`.
2. **Python stack** — the `tumult/` synthesis library and the `TumultStudio2/`
   Bokeh GUI built on top of it.

---

## 1. Python modules

**Recommended Python version: 3.11** (the reference `tumult_studio` environment
runs 3.11.13). Newer versions (3.12) also work.

All Python dependencies for `tumult/` and `TumultStudio2/` are listed in
[`requirements.txt`](requirements.txt) and installed with pip:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Launching TumultStudio2

From the repository root (the directory containing this file):

```bash
bokeh serve --show TumultStudio2
```

The `tumult/` package is imported by adding the repo root to `sys.path`, so run
Bokeh from here rather than from inside `TumultStudio2/`.

---

## 2. Compiling the Fortran code

The build uses the top-level DISPATCH `config/Makefile` (two levels up) together
with the local [`Makefile`](Makefile). The compiler is selected automatically:

- **macOS (Darwin):** `gfort10` (a GCC 10 gfortran; **not** the default `gfortran`)
- **Linux:** `ifort`

Build with:

```bash
make -j
```

This produces the executable `dispatch.x` in this directory. To clean:

```bash
make clean
```

### Running

```bash
./dispatch.x demo.nml
```

where `demo.nml` (or any of the `*.nml` namelist files here) configures the run.
Output spectra are written to the paths given in the namelist.

---