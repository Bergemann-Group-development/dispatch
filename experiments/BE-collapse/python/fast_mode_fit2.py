import matplotlib.pyplot as pl
import numpy as np
from scipy.optimize import curve_fit
from read_file import read_file
from image import image

def fun(t,a,delta,w,phi):
    return a*np.exp(-t*delta)*np.sin(w*t+phi)

def lim(y):
    pl.xlim((-1,24))
    pl.ylim(y)
    
fsz=(12,8)
pl.figure(1,figsize=fsz)
pl.show()
pl.ion()
ampl=0.01
y=(-ampl,ampl)
px=[]
for i in range(200):
    file='data/fast_mode_001_{:04d}'.format(i)
    print file
    f=read_file(file)
    sz=f.d.shape
    iz=sz[2]/4
    iy=sz[1]/2
    ix=5
    px.append(f.px[ix,iz,iz])

pl.plot(px,'-o')

nt=np.size(px)
t=np.linspace(0.0,2.0,nt)
p, pcov = curve_fit(fun, t, px)
pl.plot(fun(t,p[0],p[1],p[2],p[3]),'r')
k=1
w=p[2]; delta=p[1]
speed=p[2]/(1.0*np.pi*k)
print "        wave speed: {}".format(speed)
print "damping per period: {}".format(delta/speed)
