import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu

def yx_midplane(run='',dir='../data',iout=0,var='d'):
    '''
    Display midplane (z=0) xy-images
    '''
    pp=du.patches(run=run,dir=dir)
    pp=[p for p in pp if p.iout==iout]
    pp=[p for p in pp if p.pos[2]==0.0 and p.level<8]
    print np.size(pp),'patches in midplane'    
    pl.figure(figsize=[8,3])
    #pl.clf()
    vmin=1e9
    vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.array(p.nghost)
    u=l+np.array(p.n)
    m=(l+u)/2
    for p in pp:
        p.cache()
        p.vars()
        p.f=p.var[var]
        p.f=p.f[l[0]:u[0],l[1]:u[1],m[2]]
        e=p.extent[2]
        p.f,e=du.rotate(p.f,e)
        p.extent[2]=e
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
        #print p.id,p.level
    print vmin,vmax    
    for p in pp:
        pu.image(p.f,extent=p.extent[2],vmin=vmin,vmax=vmax)    
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.title(var+'   t={:d}'.format(np.int(p.time)))
    pl.xlabel('y')
    pl.ylabel('x')
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.gcf().set_size_inches(12, 3)
    pl.tight_layout()
    pl.show()

iout=30
run='m10r40t10'
dir='d:\\pebble_atmosphere\\data'
vars=['logd','ux','uy','uz','loge']
for var in vars:
    yx_midplane(run,dir,iout=iout,var=var)
