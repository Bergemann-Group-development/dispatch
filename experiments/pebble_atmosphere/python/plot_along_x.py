import numpy as np
reload(np)
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu
reload(pu)

''' 
    Use zoom-in in the Python graphics backend to focus on specfic
    ranges, especially when using "plot" rather than "semilogx"
'''

data='d:\\pebble_atmosphere\\data'
run='m3r27t05st'
iout=20
sgn=1

pp=du.patches(run,data,iout=iout)
pp=[p for p in pp if p.pos[1]==0.0 and p.pos[2]==0.0 and sgn*p.pos[0]>=0.0]
npp=np.size(pp)
print npp,'patches to plot'
dx0=pp[npp-1].dx[0]

pl.figure(1)
pl.show()
#pl.clf()
for p in pp:
    print p.id
    p.vars()
    l=0
    n=p.gn[0]-1
    u=l+n
    iy=(l+u)/2-1; iz=iy
    x=p.corners[0][0]+(l+np.arange(n)+0.5)*p.dx[0]
    x=sgn*x/dx0
    pg =(p.data[l:u,iy  ,iz  ,4]\
        +p.data[l:u,iy+1,iz  ,4]\
        +p.data[l:u,iy  ,iz+1,4]\
        +p.data[l:u,iy+1,iz+1,4])/4
    pg=pg*(p.gamma-1)
    d  =(p.data[l:u,iy  ,iz  ,0]\
        +p.data[l:u,iy+1,iz  ,0]\
        +p.data[l:u,iy  ,iz+1,0]\
        +p.data[l:u,iy+1,iz+1,0])/4
    px = p.p1[l:u,iy,iz]
    ux =(p.u1[l:u,iy  ,iz  ]\
        +p.u1[l:u,iy+1,iz  ]\
        +p.u1[l:u,iy  ,iz+1]\
        +p.u1[l:u,iy+1,iz+1])/4
    s  =(p.s[l:u,iy  ,iz  ]\
        +p.s[l:u,iy+1,iz  ]\
        +p.s[l:u,iy  ,iz+1]\
        +p.s[l:u,iy+1,iz+1])/4
    i0=0
    for i in range(n):
        if sgn*x[i]<=0.0:
            i0=i+1
    if sgn>0:
        x=x[i0:n]
        s=s[i0:n]
        d=d[i0:n]
        ux=ux[i0:n]
        px=px[i0:n]
    else:
        x=x[0:i0-1]
        s=s[0:i0-1]
        d=d[0:i0-1]
        ux=ux[0:i0-1]
        px=px[0:i0-1]
    pl.subplot(2,2,1)
    pl.semilogx(x,px,'-+')
    pl.title('{}  iout={}  (:,{},{})'.format('px',iout,iy,iz))
    pl.subplot(2,2,2)
    pl.semilogx(x,ux,'-+')
    pl.title('{}  iout={}  (:,{},{})'.format('ux',iout,iy,iz))
    pl.subplot(2,2,3)
    #pl.plot(x,s ,'-+')
    pl.semilogx(x,s ,'-+')
    pl.title('{}  iout={}  (:,{},{})'.format('s',iout,iy,iz))
    pl.subplot(2,2,4)
    pl.loglog(x,d ,'-+')
    pl.title('{}  iout={}  (:,{},{})'.format('d',iout,iy,iz))
pl.draw()
