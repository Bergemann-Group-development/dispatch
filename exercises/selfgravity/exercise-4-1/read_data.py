'''
   Read the density and potential from the data directory
'''
import numpy as np
import matplotlib.pyplot as pl

f = open('data/poisson_001_0000.txt','r')
# line 1: time pos(1:3) size(1:3)
line = f.readline()
column = line.split()
time = float(column[0])
pos = [float(column[1]),float(column[2]),float(column[3])]
size = [float(column[4]),float(column[5]),float(column[6])]
print "           time:",time
print "       position:", pos
print "           size:", size
# line 2: m(1:3) n(1:3) li(1:3) ui(1:3)
line = f.readline()
column = line.split()
m = [int(column[0]),int(column[1]),int(column[2])]
n = [int(column[3]),int(column[4]),int(column[5])]
li = [int(column[6]),int(column[7]),int(column[8])]
ui = [int(column[9]),int(column[10]),int(column[11])]
print "     full patch:", m
print "  partial patch:", n
print "  lower indices:", li
print "  upper indices:", ui
f.close 

f = open('data/poisson_001_0000.dat','rb')
shape = (n[0],n[1],n[2])
rho = np.fromfile(file=f, dtype=np.float32, count=np.prod(shape)).reshape(shape)
print "  density range:", [rho.min(),rho.max()]
pl.subplot(1,2,1)
pl.imshow(rho[:,:,n[2]/2],origin='lower',interpolation='nearest')
pl.title('density')
phi = np.fromfile(file=f, dtype=np.float32, count=np.prod(shape)).reshape(shape)
print "potential range:", [phi.min(),phi.max()]
pl.subplot(1,2,2)
pl.title('potential')
pl.imshow(phi[:,:,n[2]/2],origin='lower',interpolation='nearest')
pl.show()
