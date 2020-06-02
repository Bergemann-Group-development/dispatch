import os
import sys
import numpy as np
import matplotlib.pyplot as pl
import dispatch as dis
import dispatch.select as dse
import dispatch.graphics as dgr
from sink_data import sink_data

run='128'
if len(sys.argv) > 1:
    run=sys.argv[1]
data='../data'
dir=data+'/'+run+'/figures'
if not os.path.exists(dir): os.mkdir(dir)
v=sink_data(run)

sn=dis.snapshot(0,run)
m=v.m*sn.scaling.m/sn.cgs.m_sun
t=v.t*sn.scaling.t/sn.cgs.myr
pl.figure(); pl.plot(t,m)
pl.xlabel('t [Myr]'); pl.ylabel('m [solar]')
pl.savefig(dir+'/tm.png')
pl.close()

x=v.x*sn.scaling.l/sn.cgs.pc
y=v.y*sn.scaling.l/sn.cgs.pc
pl.figure(); pl.plot(x,y)
for io in range(9,51):
    t0=io*0.01
    w=np.where(v.t > t0)[0][0]
    pl.plot(x[w],y[w],'+',color='b')
pl.xlabel('x [pc]'); pl.ylabel('y [pc]')
pl.savefig(dir+'/xy.png')
pl.close()
