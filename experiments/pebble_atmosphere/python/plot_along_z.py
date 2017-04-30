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
run='m3r27t05'
iout=5
sgn=1

pp=du.patches(run,data,iout=iout)
pp=[p for p in pp if p.pos[0]==0.0 and p.pos[1]==0.0 and sgn*p.pos[2]>=0.0]
npp=np.size(pp)
print npp,'patches to plot'
dz0=pp[npp-1].dx[2]

pl.figure(1)
pl.show()
#pl.clf()
for p in pp:
    p.vars()
    l=0
    n=p.gn[2]
    u=l+n
    iy=(l+u+1)/2; ix=iy
    z=p.corners[0][2]+(l+np.arange(n)+0.5)*p.dx[2]
    z=sgn*z #/dz0
    if p.id==60:
        p60=p
        z60=z
    elif p.id==141:
        p141=p
        z141=z
    elif p.id==78:
        p78=p
        z78=z
    elif p.id==159:
        p159=p
        z159=z
    #print p.id,z.max()
    pg =(p.data[ix  ,iy  ,l:u,4]\
        +p.data[ix+1,iy  ,l:u,4]\
        +p.data[ix  ,iy+1,l:u,4]\
        +p.data[ix+1,iy+1,l:u,4])/4
    pg=pg*(p.gamma-1)
    d  =(p.data[ix  ,iy  ,l:u,0]\
        +p.data[ix+1,iy  ,l:u,0]\
        +p.data[ix  ,iy+1,l:u,0]\
        +p.data[ix+1,iy+1,l:u,0])/4
    pz = p.p3[ix,iy,l:u]
    uz =(p.u3[ix  ,iy  ,l:u]\
        +p.u3[ix+1,iy  ,l:u]\
        +p.u3[ix  ,iy+1,l:u]\
        +p.u3[ix+1,iy+1,l:u])/4
    s  =(p.s[ix  ,iy  ,l:u]\
        +p.s[ix+1,iy  ,l:u]\
        +p.s[ix  ,iy+1,l:u]\
        +p.s[ix+1,iy+1,l:u])/4
    i0=0
    for i in range(n):
        if sgn*z[i]<=0.0:
            i0=i+1
    if sgn>0:
        z=z[i0:n]
        s=s[i0:n]
        d=d[i0:n]
        uz=uz[i0:n]
        pz=pz[i0:n]
    else:
        z=z[0:i0-1]
        s=s[0:i0-1]
        d=d[0:i0-1]
        uz=uz[0:i0-1]
        pz=pz[0:i0-1]
    pl.subplot(2,2,1)
    zs=z-sgn*0.499*p.dx[2]
    if p.id==60:
        zs60=zs
    elif p.id==141:
        zs141=zs
    print p.id,i0,np.size(z),zs[0:3],'  ',zs[n-3:n]
    pl.semilogx(zs,pz,'-+')
    pl.title('{}  t={}  ({},{},:)'.format('pz',p.time,ix,iy))
    pl.subplot(2,2,2)
    pl.semilogx(z,uz,'-+')
    pl.title('{}  t={}  ({},{},:)'.format('uz',p.time,ix,iy))
    pl.subplot(2,2,3)
    #pl.plot(z,s ,'-+')
    pl.semilogx(z,s ,'-+')
    pl.title('{}  t={}  ({},{},:)'.format('s',p.time,ix,iy))
    pl.subplot(2,2,4)
    pl.loglog(z,d ,'-+')
    pl.title('{}  t={}  ({},{},:)'.format('d',p.time,ix,iy))
pl.draw()
