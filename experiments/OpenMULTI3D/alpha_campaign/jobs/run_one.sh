#!/bin/bash
# run_one.sh NML THREADS -- one dispatch.x run through common.sh's run(), for xargs.
source /mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign/jobs/common.sh
run "$1" "${2:-1}"
