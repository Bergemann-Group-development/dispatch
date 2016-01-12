'''
   Read particle data from an "unformatted" Fortran file.  A record
   in such a file has a specific structure, which the FortranFile
   package can interprete. The f.readInts() and f.readReals() methods
   read as many ints and reals as possible, respectively.

   As is, this procedure reads in the data for all particles (which
   could be one or many), and plots the projection onto the x-y plane
   of the first particle.  You can modify the procedure, e.g. to plot
   more orbits, or do other computations on the data, by adding source
   text at the bottom.
'''
import numpy as np
import matplotlib.pyplot as pyplot
from fortranfile import FortranFile

all = []

with FortranFile("data/kdk.dat", "r") as f:
    # Each loop in the while tries to read first an integer
    # giving the number of particles that follow, and then 
    # reads the particle data, which for each particls is
    # time,x,y,z,vx,vy,vz,mass.  These are collected into the
    # one-dimensional array "all", as long as there is more
    # to read, which could be for a number of time steps
    j = 0
    while True:
        try:
            n = f.readInts()
            n = n[0]
        except:
            break
        for i in range(n):
            p = f.readReals()
            print j, p[0]
            all.append(p)
            j = j+1

# After finished reading, the data is reshaped in to an array
# with dimensions (ntime,nparticles,8), where index 0 is time,
# indices 1-3 # are {x,y,z}, indices 4-6 are velocities, and
# index 7 is mass
all = np.concatenate(all)
s = all.shape
all.shape = (s[0]/n/8,n,8)

# plot the x and y coordinates of the first particle
pyplot.plot(all[:,0,1],all[:,0,2],'o')
pyplot.show()

# Additional processing may be added here:
