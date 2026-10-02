#!/bin/bash
# Batch 2 (Saturday): the paired pilot on seven Stagger models with K snapshots each
# (default 2): 1D LTE/NLTE on the MARCS twins and <3D> averages, 3D LTE, 3D NLTE Mg,
# and (separately, on request) 3D NLTE Fe.  Walltimes from the phase-0 timings with a
# factor ~2 margin; see reports/2026-10-02_phase0_recon.md.
set -e
K=${1:-2}
CAMPAIGN=/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign
cd "$CAMPAIGN/work"
python3 ../jobs/make_nml.py pilot "$K"
cd "$CAMPAIGN"
# 42 1D runs, 8 threads each, 4 at a time: ~5 min
sbatch -t 00:30:00 -J p1_1d jobs/run_tasks.sbatch p1/tasks_1d.txt 8
# 3D LTE, 7K runs at 32 threads, ~15 min each at R=2e5
sbatch --array=0-1 -t 06:00:00 -J p1_3dlte jobs/run_tasks.sbatch p1/tasks_3d_lte.txt 32
# 3D NLTE Mg with its NLTE and LTE spectra, 7K runs at 32 threads, ~1.5 h each (bench:
# departures 64 min on the solar snapshot, 3 GiB)
sbatch --array=0-1 -t 24:00:00 -J p1_3dmg jobs/run_tasks.sbatch p1/tasks_3d_mg.txt 32
# 3D NLTE Fe (~80 GB, hours per snapshot): submit by hand once the Mg runs are checked,
#   sbatch --array=0-6 -t 48:00:00 -J p1_3dfe jobs/run_tasks.sbatch p1/tasks_3d_fe.txt 32
