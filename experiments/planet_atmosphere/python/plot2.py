from read_file import read_file
import matplotlib.pyplot as pl
import numpy as np

buf=np.zeros((30,19),dtype=np.float32)
for i in range(19):
    file='data/planet_atmosphere_068_{:04d}'.format(i)
    f=read_file(file)
    buf[:,i]=f.px[:,15,15]
