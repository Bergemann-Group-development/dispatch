import numpy as np
import matplotlib.pyplot as pl
import fortranfile as ftf

file='fort.22'
ee=[]
with ftf.FortranFile(file,'rb') as f:
    n=0
    while True:
        try:
            e=f.readReals(prec='d')
            ee.append(e)
        except:
            break
        n=n+1 
print np.shape(ee)
ee=np.array(ee)

pl.figure(1)
pl.clf()
pl.plot(ee[:,0],1-ee[:,4]/ee[0,4])
pl.show()

'''
file='fort.32'
with ftf.FortranFile(file,'rb') as f:
    i=f.readInts()
    n_period=i[0]
'''