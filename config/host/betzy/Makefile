COMPILER = ifort

define HOSTINFO

  Best to use this modules: (last update 30.04.2021)

    module restore
    module load iompi/2020a

  SLURM script example: `betzy.run`

    #!/bin/bash

    #SBATCH --job-name=dispatch
    #SBATCH --time=0-12:00:00
    #SBATCH --account=nn2834k
    #SBATCH --comment="Dispatch/Bifrost Stellar Atmosphere"

    ## for pure mpi, one task per core:
    ## set nodes only. In example: -np 16384 = 128x128

    #SBATCH --nodes=4

    ## Set one mpi process per socket/cpu
    #SBATCH --ntasks-per-node=2

	 ## Set 64 threads per MPI rank     
    #SBATCH --cores-per-socket=64

    #SBATCH --output="%x_%u_%j.out"
    #SBATCH --error="%x_%u_%j.err"

    ## You can turn on notification:
    ##SBATCH --mail-type=ALL
    ##SBATCH -M my.email@astro.uio.no

    ## Type of job sorted by priority
    ##SBATCH --qos=devel   # max nodes 4, max time 30 min
    ##SBATCH --qos=preproc # max nodes 1, max time 1 day
    ##SBATCH --qos=short   # max nodes 10, max time 2 hours
    ##SBATCH --qos=normal  # max time one week

    module -q restore
    module load iompi/2020a

    ## Recommended safety settings:
    set -o errexit # Make bash exit on any error
    set -o nounset # Treat unset variables as errors

    cd $$SLURM_SUBMIT_DIR

    now=`date '+%Y-%m-%d_%H:%M:%S'`
    log_file="$$(basename $$(pwd))_$$now.log"

    # Optional setup to increase stack 
    ulimit -s unlimited
	 export OMP_STACKSIZE=512M

    cat $$0 | tee $$log_file
    time srun --mpi=pmix ./dispatch.x test.nml | tee -a $$log_file

    exit 0

  Submit with:

    sbatch betzy.run

  Check progress with:

    squeue -u username

  Cancel with:

    scancel job_id

endef
