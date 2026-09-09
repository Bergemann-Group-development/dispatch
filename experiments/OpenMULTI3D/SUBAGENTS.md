# Working in this repository as a subagent

Read `AGENTS.md` before starting. Your checkout is already a clean worktree on a
task-specific branch based on the shared checkpoint named in your assignment.

Commit the completed task to your branch in one focused commit. Report the commit hash, the purpose of the changed files, verification commands with
their actual results, and any remaining limitation. Do not integrate or remove the
branch or worktree yourself.

## Repository-specific operations

`SRC/*.i90` files are compiler preprocessor output. Edit the corresponding `.f90`
source.

Run the executable as `./dispatch.x <namelist>` from the project directory. The
namelist argument must be a short relative path. A long absolute path overflows a
fixed-length buffer in the parent `io_mod` and aborts with
`Fortran runtime error: End of record` during startup. That abort is an artifact of
how the program was invoked. `datadir` inside the namelist may be a long absolute
path; only the argument is constrained. Remove scratch namelists from the project
directory when finished.

The Python environment for tests and analysis is the `tumult_studio` conda
environment. Establish the test count before a change when comparing suite results so
that any difference can be attributed correctly.

## Verification

Useful techniques in this repository:

- A run writes its outputs under `datadir/<run>/save/`. Byte-comparing every `.bin`
  there against a baseline run is a strong regression check for a change that should
  not alter numerics. Produce the baseline from an unmodified build rather than
  assuming what it would contain.
- `save/flux_meta.txt` lists one row per output patch with its point count, which is how
  to confirm that a change to patch decomposition took effect.
- A control run that exercises the existing path separates a disabled feature from a
  failure in the underlying machinery.
- Read spectra through `tumult.read()` and `crop` rather than parsing saved binaries.
  Comparing a raw wavelength array directly is reasonable; say so when doing it.
- Search the run log for warnings the change can emit. A guard that fires when it
  should not may be invisible in a run that otherwise completes normally.
