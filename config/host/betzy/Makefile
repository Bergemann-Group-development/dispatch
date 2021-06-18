COMPILER=intel

define HOSTINFO

  Best to use intel (hiden) modules: (last update 10.05.2021)

    module restore
    module --ignore-cache load intel/2020b

  SLURM script example: `betzy.run`

    #!/bin/bash

    #SBATCH --job-name=dispatch
    #SBATCH --time=0-12:00:00
    #SBATCH --account=nn2834k
    #SBATCH --comment="Dispatch/Bifrost Stellar Atmosphere"

    ## for best performance there is a fixed number of 8 MPI process
    ## per node, thus set nodes only. In example 4*8=32 == -np 32

    #SBATCH --nodes=4

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
    module --ignore-cache load intel/2020b

    ## Recommended safety settings:
    set -o errexit # Make bash exit on any error
    set -o nounset # Treat unset variables as errors

    cd $$SLURM_SUBMIT_DIR

    now=`date '+%Y-%m-%d_%H:%M:%S'`
    log_file="$$(basename $$(pwd))_$$now.log"

    ## Optional setup to increase stack
    ulimit -s unlimited
	  export OMP_STACKSIZE=256M

    ## MPI & OMP processes pinning
    export I_MPI_PIN_DOMAIN=numa
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread

    mpirun -genv OMP_NUM_THREADS=8 -ppn 8 ./dispatch.x test.nml | tee $$log_file

    cat $$0 | tee $$log_file
    exit 0

  Submit with:

    sbatch betzy.run

  Check progress with:

    squeue -u username

  Cancel with:

    scancel job_id

endef
