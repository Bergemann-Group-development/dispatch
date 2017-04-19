import numpy as np
import matplotlib.pyplot as pl
import dispatch_utils as du
reload(du)
import plot_utils as pu

def yz_midplane(run='',dir='../data',iout=0,var='d'):
    '''
    Display midplane (x=0) yz-images
    '''
    pp=du.patches(run=run,dir=dir)
    pp=[p for p in pp if p.iout==iout]
    pp=[p for p in pp if p.pos[0]==0.0]
    print np.size(pp),'patches in midplane' 
    pl.figure(figsize=[8,3])
    pl.clf()
    vmin=1e9
    vmax=-vmin
    xlim=[vmin,vmax]
    ylim=[vmin,vmax]
    l=np.array(p.nghost)
    u=l+np.array(p.n)-1
    m=(l+u)/2
    for p in pp:
        p=du.cache_patch(p)
        if var=='d':
            p.f= p.d[m[0],l[1]:u[1],l[2]:u[2]]
            p.f=np.log10(p.f)
        elif var=='ux':
            p.f=p.ux[m[0],l[1]:u[1],l[2]:u[2]]
        elif var=='uy':
            p.f=p.uy[m[0],l[1]:u[1],l[2]:u[2]]
        elif var=='uz':
            p.f=p.uz[m[0],l[1]:u[1],l[2]:u[2]]            
        elif var=='e':
            p.f= p.e[m[0],l[1]:u[1],l[2]:u[2]]            
            p.f=np.log10(p.f)
        e=p.extent[0]
        vmin=min(vmin,p.f.min())
        vmax=max(vmax,p.f.max())
        xlim=[min(xlim[0],e[0]),max(xlim[1],e[1])]
        ylim=[min(ylim[0],e[2]),max(ylim[1],e[3])]
    print vmin,vmax   
    for p in pp:
        pu.image(p.f,extent=p.extent[0],vmin=vmin,vmax=vmax)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.title(var+'   t={:d}'.format(np.int(p.time)))
    pl.xlabel('y')
    pl.ylabel('z')
    pl.colorbar(shrink=0.8,aspect=15,fraction=0.04,pad=0.02)
    pl.gcf().set_size_inches(12, 3)
    pl.tight_layout()
    pl.show()

iout=2
run='m10r40t10'
dir='d:\\pebble_atmosphere\\data'
vars=['d','ux','uy','uz','e']
for var in vars:
    yz_midplane(run,dir,iout=iout,var=var)
