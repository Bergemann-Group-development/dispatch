'''
   Demonstrating the consquences of the different index order in Fortran and Python
'''
import numpy as np

# Step-by-step: reading into a 1-D array, reshaping into 3-D, and transposing into Fortran order
f=open('transpose.dat','rb')
shape=(2,3,4)
a=np.fromfile(file=f, dtype=np.float32, count=np.prod(shape))
a=a.reshape(shape)
a=a.transpose(2,1,0)
f.close
    
for iz in range(0,2):
    print a[:,:,iz]
    print ''
    
# Compact form: direct application of reshape and transpose
f=open('transpose.dat','rb')
shape=(2,3,4)
a=np.fromfile(file=f, dtype=np.float32, count=np.prod(shape)).reshape(shape).transpose(2,1,0)
f.close
    
for iz in range(0,2):
    print a[:,:,iz]
    print ''
    