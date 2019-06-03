
jobcontrol.csh
==============

The job script hierarchy consists of SUBMIT scripts,
JOB scripts, and CANNED scripts.

The SUBMIT scipts contain scheduler control commands,
but other otherwise static, since they are only read at
submission.

The JOB scripts contain whatever the user wants, and may
be changed during the runs -- e.g. before restarting a
process via jobcontrol.csh.

The CANNED scripts are helpful with for example archiving
copies of executables, input files and log files.

Submit scripts:
---------------
::

        long.sh
        short.sh
        ...

These contains control syntax such as #PBS or #SBATCH.
Since changes after submission have no effect, they call
other scripts where changes do have effect.  The call is
via the jobcontrol.csh script, with submit script syntax
such as

source $HOME/codes/dispatch/config/startup.sh
jobcontrol.csh job.sh

The jobcontrol.csh script runs the script given as argument
repeatedly, as long as there is a restart.flag file present
in the run directory.

Job scripts:
------------
::

        long.csh
        short.csh

These can be anything, including calls to standard scripts for
local systems, with syntax such as

ivy.sh input.nml

Canned scripts:
---------------
::

        config/host/steno/ivy.sh

Since we often want to use standard ways of handling executables,
input files, and log file, we keep a number of canned scripts,
such as ivy.sh, which places a copy of the excutable, the input
file, and the log file, in the data/run/ directory.

The script seach path is setup by the config/startup.sh script,
with directories seached in this order::

        ./
        config/host/$HOST/
        utilities/scripts/$SCHEDULER/
        utilities/scripts/
        $PATH

This meane one can override the default scripts by placing scripts
with the same name (e.g. ivy.sh) in the run directory.

SCHEDULER is set by a file config/host/$HOST/env.{csh,sh}, if it
exists, or may be set in the environment.  Possible values are
"pbs", "slurm", etc.

The utilities/scripts/$SCHEDULER/ directories contain scheduler-i
specific scripts, used e.g. by jobcontrol.csh.  Some examples are::

    jobid           # return the id of the current job
    jobkill id      # cancel a job
    jobmail         # send a mail to $USER

.. toctree::
   :maxdepth: 3

