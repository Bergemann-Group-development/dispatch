"""
    Python script to extract parameters from all data/run*/log_* files, where "run" is the
    argument to the script, to be executed with for example "python python/params.py 1024".

    With detailed comments, so it may be used as a an example on how to handle text in Python.
"""
import os
import sys
import numpy as np

# pick up the command line argument (the 0-th argument is the name of the script)
dir=sys.argv[1]

# Read the output from a unix command
cmd = "egrep 'snaps|PK|T_TURN|AMPL' data/"+dir+"*/log_*"
print 'extracting info from command: '+cmd
with os.popen(cmd) as fd:
    flag=False
    nt=0
    for line in fd:
        # Split each line into columns by white space
        cols = line.split()

        # Detect lines that give the PK parameter
        if cols[1]=='PK':
            if flag:
                # If we have collected info from a previous run, print it out
                pr()
            # Convert the column 3 text to float, after stripping a comma
            pk = np.float(cols[3].strip(','))

        # Detect lines that give the AMPL_TURB parameter
        if cols[1]=='AMPL_TURB':
            ampl = np.float(cols[3].strip(','))

        # Detect lines that give the T_TURN parameter
        if cols[1]=='T_TURN':
            t_turn = np.float(cols[3].strip(','))
            words=cols[0].split('/')
            run=words[1]
            log=words[2].strip(':')

        # Detect lines that give snapshot number and time.  Since these lines follow after the
        # lines that give parameters, we raise a flag, and delay printing until a line with a
        # parameter from a new run is found, or the input ends
        if cols[1]=='snapshot':
            flag=True

            # pick up the snapshot time
            t = np.float(cols[4])

            # split the column data/dir/log_xxx in /, and pick up the snapshot number
            words=cols[2].split('/')
            nt = np.int(words[2])

        # Now that all information has been gathered, a function can be defined that prints everything, when called (later)
        def pr():
            print '{}   {}   ampl: {:.1f}    t_turn: {:.4f}    pk: {:.2f}   nt: {:3d}   t: {:.3f}'.format(run,log,ampl,t_turn,pk,nt,t)
            flag=False

# Printo out the info from the last run
pr()
