# STEPS=vmic for queue.sbatch: the vmic-grid 1D LTE partners of the 3D LTE snapshots.
# The list is regenerated when the submission starts (make_nml.py vmic_grid), so that a
# job queued behind the avg3d queue sees the 3D LTE runs and <3D> averages finished by then.
#   sbatch -c 1 --mem=3G --export=ALL,STEPS=vmic jobs/queue.sbatch b5/tasks_vmic.txt
prepare () { python3 $CAMPAIGN/jobs/make_nml.py vmic_grid; }
