import numpy as np
import matplotlib.pyplot as pyplot
from fortranfile import FortranFile

all = []

a = np.arange(2, dtype='int8')
with FortranFile("data/particles.dat", "r") as f:
  while True:
    try:
      a = f.readInts()
    except:
      break
    n = a[0]
    p = f.readReals()
    all.append(p)
  
all = np.concatenate(all)
s = all.shape
all.shape = (s[0]/8,8)
pyplot.plot(all[:,1],all[:,2],'o')
pyplot.show()