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

data='../data'
run='3me27r'
iout=0
sgn=-1
#pl.interactive(True)
pp=du.patches(iout,run,data)
pp=[p for p in pp if p.pos[0]==0.0 and p.pos[2]==0.0 and sgn*p.pos[1]>=0.0]
npp=np.size(pp)
print npp,'patches to plot'
dy0=pp[npp-1].dx[2]

pl.figure(1)
#pl.clf()
for p in pp:
    p.vars()
    l=0
    n=p.gn[1]
    u=l+n
    iz=12; ix=iz-1
    y=p.corners[0][1]+(l+np.arange(n)+0.5)*p.dx[1]
    y=sgn*y #/dz0
    if p.id==60:
        p60=p
        y60=y
    elif p.id==141:
        p141=p
        y141=y
    elif p.id==78:
        p78=p
        y78=y
    elif p.id==159:
        p159=p
        y159=y
    #print p.id,z.max()
    pg =(p.data[ix  ,l:u,iz  ,4]\
        +p.data[ix+1,l:u,iz  ,4]\
        +p.data[ix  ,l:u,iz+1,4]\
        +p.data[ix+1,l:u,iz+1,4])/4
    pg=pg*(p.gamma-1)
    d  =(p.data[ix  ,l:u,iz  ,0]\
        +p.data[ix+1,l:u,iz  ,0]\
        +p.data[ix  ,l:u,iz+1,0]\
        +p.data[ix+1,l:u,iz+1,0])/4
    py = p.p3[ix,l:u,iz]
    if p.id ==39:
		print py
    uy =(p.u3[ix  ,l:u,iz]\
        +p.u3[ix+1,l:u,iz]\
        +p.u3[ix  ,l:u,iz+1]\
        +p.u3[ix+1,l:u,iz+1])/4
    s  =(p.s[ix  ,l:u,iz  ]\
        +p.s[ix+1,l:u,iz  ]\
        +p.s[ix  ,l:u,iz+1]\
        +p.s[ix+1,l:u,iz+1])/4
    i0=0
    for i in range(n):
        if sgn*y[i]<=0.0:
            i0=i+1
    if sgn>0:
        y=y[i0:n]
        s=s[i0:n]
        d=d[i0:n]
        uy=uy[i0:n]
        py=py[i0:n]
    else:
        y=y[0:i0-1]
        s=s[0:i0-1]
        d=d[0:i0-1]
        uy=uy[0:i0-1]
        py=py[0:i0-1]
    pl.subplot(2,2,1)
    ys=y-sgn*0.499*p.dx[1]
    if p.id==60:
        ys60=ys
    elif p.id==141:
        ys141=ys
    print p.id,i0,np.size(y),ys[0:3],'  ',ys[n-3:n]
    pl.semilogx(ys,py,'-+')
    pl.title('{}  t={}  ({},:,{})'.format('pz',p.time,ix,iz))
    pl.subplot(2,2,2)
    pl.semilogx(y,uy,'-+')
    pl.title('{}  t={}  ({},:,{})'.format('uz',p.time,ix,iz))
    pl.subplot(2,2,3)
    #pl.plot(z,s ,'-+')
    pl.semilogx(y,s ,'-+')
    pl.title('{}  t={}  ({},:,{})'.format('s',p.time,ix,iz))
    pl.subplot(2,2,4)
    pl.loglog(y,d ,'-+')
    pl.title('{}  t={}  ({},:,{})'.format('d',p.time,ix,iz))
pl.draw()
pl.show()
