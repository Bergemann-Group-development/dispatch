# Alpha campaign: local operations notes

Facts about the MPIA cluster, the campaign data and the job scripts that do not belong in
the user manual. Sources: the phase 0-5 reports in `reports/`. Statements marked "unconfirmed" come
from the wishlist only.

## Layout

- `$C = /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign` on the cluster.
  `$C/bin/<hash>` holds frozen builds, `$C/work` the namelists and run folders, `$C/logs`
  the job logs, per-run timings and submission records (`logs/submitted_*.txt`).
- `jobs/common.sh` sets `CAMPAIGN`, `WORK` and `BIN` (default `bin/current`, overridable
  with `--export=ALL,BIN=$C/bin/<hash>`) and refuses to run outside nodes 81-90. It prints
  one `binary:` line at the start of the job output and as the first line of every run log
  (`logs/<run>.log`): the resolved directory, the commit from the `REVISION` file that
  `jobs/build.sbatch` writes into each frozen build (older builds: the directory name) and
  the md5 of `dispatch.x`. Check that line before trusting a result; a job that does not
  export `BIN` runs whatever `bin/current` points to when it starts.
- Frozen builds are named by the commit hash and report it in the `DISPATCH code framework`
  banner and in `job.nml` (`DISPATCH_HASH`). The Makefile rebuilds `io_mod.o` and
  `job_record_mod.o` when HEAD moves; before that rule, an incremental build in `$C/src`
  kept the hash of the first commit it had compiled (`bin/4890f9fde` printed `f63e2fc9`).
  Builds before this change can carry a stale banner hash; the directory name and
  `REVISION` are the reliable record.
- `jobs/queue.sbatch LIST [THREADS]` runs task lists (see its header). Claims live in
  `logs/claims_<list>`; the older runners used `logs/claims`, `claims_vmic`, `claims_dep4`
  and `claims_dep4_b6`, so pass `CLAIMS=` to continue one of their lists. Pipelines per
  line (departures then SED, a list regenerated at start) are step files in `jobs/steps/`.
- Metal-poor absmet tables are in `$C/work/absmet/`; the solar table is
  `input/absmet_files/absmet`.

## Cluster quirks

- **Login shell.** Jobs start with `#!/bin/bash -l` so that the Intel runtime is on the path.
- **node83.** `/opt/intel/oneapi` is a stale NFS handle there, so ifort-built binaries fail
  with `libiomp5.so: cannot open shared object file`. `jobs/queue.sbatch` puts
  `$C/lib/libiomp5.so` (copied from the login node's oneAPI 2024.0) on `LD_LIBRARY_PATH`.
  Python jobs run unchanged. The mount is still open with the administrators.
- **Partial nodes.** Other users hold cores and memory on some nodes (node83 had four free
  cores, node88 one free core and about 5 GB, node90 no free memory during phase 3). Size a
  job with `-c` and `--mem` against `scontrol show node` first.
- **Stack size.** `KMP_STACKSIZE=256M` is set in `common.sh`; HOP from Python needs
  `ulimit -s unlimited`.
- **Allowed nodes.** Campaign jobs name `--nodelist=mpia-slurm-node[81-90]`. A job without a
  node list landed eight tasks on node47.
- **Node rule** (replaces the earlier weekday rule of `--nice` and at most four whole
  nodes). On nodes 81-90 (87 is down, 89 does not exist):
  1. Always leave at least one usable node completely idle, counting other users' jobs.
     A node that is down or drained does not count as idle.
  2. Claim whole nodes (`-c` equal to the node's cores, memory to match), never half a node.
  3. A node that our jobs already partly occupy may be filled up.
  Before submitting, count the idle nodes: `nodes-free` (skill `mpia-cluster`) lists each
  node with its cores held by us and by others, and says whether one more whole node may be
  claimed. With a single idle node left, fill a partly used node of ours or wait.

## 4MOST data

- Spectra: `~/Documents/Spectra/4MOST/normalised_rv_corrected_19june2026`. HRS files are
  text tables (wavelength, normalised flux, error) in three arms of 0.05 Å sampling,
  3926-4355, 5160-5730 and 6100-6790 Å (air); LRS files cover 3700-9443 Å at 0.25 Å.
  The `stack1/3/4` subfolders are undocumented.
- Sets: 0201 (675 spectra, 85 stars, HRS, mostly giants, plausibly the low-metallicity
  programme), 0301 (84, HRS giants), 0401 (463, HRS dwarfs and subgiants), 0408 (1875, HRS
  turn-off and subgiants), 0409 (photometry only, red giants with asteroseismic data), s3
  (1187, LRS).
- No pipeline labels are included. The file names say `sapp_input_spectrum_from_fgkpayne_*`;
  the labels are probably held by the 4MOST consortium pipeline.
- Photometry CSVs carry Gaia DR3, 2MASS, distances and some radial velocities. Match IDs as
  strings.
- Related material on beegfs: `storm/4MOST/` (Payne training grids and
  `victor_predictions_sep2025/` labels for benchmark-star simulations),
  `storm/payne/` (Gaia-ESO setups), `storm/data/SAPP_DR4pub` (Gaia-ESO DR4 spectra),
  `users/shared-storage/gent/4MOST_OpR2_spec_analysis/` (simulated OpR2 data and the SAPP
  code). None of these contain labels for the 2026 spectra.
- From the wishlist, unconfirmed here: the FITS columns are big-endian and `mg_gspphot` is an absolute magnitude.

## Code history worth knowing

- `&atmos_params prune_top`/`prune_tol` lower the top of a 3D box in the Fortran reader by
  an electron-pressure criterion (removed in 22e913959, restored by its revert). `nz`
  resamples the remaining box. The campaign's depth crop at a fixed <log tau_500> is the
  data-side step `jobs/crop3d.py`.
- Two native HOP tests failed before the campaign: a stale `hop_opacity` golden and no host
  caller for `hop_holtsmark`.

## Open items

See `reports/2026-10-05_agent_wishlist.md` for the ranked list of tooling and documentation
work.
