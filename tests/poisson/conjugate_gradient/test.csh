#!/bin/csh

set dir = performance
mkdir -p $dir

set host = `hostname`
set log = $dir/$host.log
echo "log file: $log"

hostname | tee $log
grep 'model name' /proc/cpuinfo | head -1 | tee -a $log
foreach m ( 'cg' 'sor' )
foreach o ( 0 1 2 3)
foreach t ( f t )
./test.x << EOF | grep TIME | tee -a $log
&cg_params
 dims=3*48
 move=6
 tolerance=1e-4
 floor=1e-2
 order=$o
 width=0.1
 method='$m'
 red_black=$t
 precondition=$t
 max_iter=900
/
EOF
end
end
end
