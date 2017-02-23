import numpy as np
import matplotlib as mp
from dispatch_data import *
from dispatch_grid import *
from derived_variables import *
import itertools
import yt as yt
import scaling as sc

ptimes=[0.002]
t_eps=0.0001
ThisPatch = RTPatch
verbose=False
skipped=[]
orbit_radius_au = 1.0
omega = [0.0,0.0,1.0]
pl_pos= [0.0,1.0,0.0]
v_headwind = 0.0   # the headwind the planet feels
v_accretion = 0.0  # accretion rate
density_cgs=3e-12  # external density
pl_rad=sc.scaling().l
for i in range(0,2):
    omega[i]=omega[i]*2.*np.pi*sc.scaling().t/sc.yr
    pl_pos[i] = pl_pos[i]*orbit_radius_au*sc.au/sc.scaling().l
v_headwind=v_headwind/sc.scaling().u
v_accretion=v_accretion/sc.scaling().u
rho_gas = density_cgs/sc.scaling().d

def get_velocity(n, rho, p):
    vel = np.zeros([n[0],n[1],n[2]])
    for ii in range(0,n[0]):
        for ij in range(0,n[1]):
            for ik in range(0,n[2]):
                if p[ii,ij,ik] !=0.:
                    vel[ii,ij,ik]=p[ii,ij,ik]/rho[ii,ij,ik]
    return vel
def get_data(time,t_eps,datadir,prefix='',verbose=False):
    files = fnmatch.filter(os.listdir(datadir),prefix+"*.txt")
    patches=[]
    times=[]
    for f in files:
        ptemp = ThisPatch(datadir+f.strip('.txt'),verbose=verbose)
        if(abs(ptemp.time-time)/(time+1.0e-99) <= t_eps):
            patches.append(ptemp)
        for p in patches:
            times.append(p.time)
    print len(patches),"patches collected for time", time
    return patches

def velocities(patches):
    for p in patches:
        
        if p.task_type == 'immersed_boundary_t':
                circle = mplpatches.Circle((p.pos[0],p.pos[1]),radius=pradius,color='k')
                cax.add_patch(circle)
                print "Skipping ",p.id," (immersed boundary)"
                continue
        
    return patches

def make_yt_dict(patches):
    gd=[]
    n=patches[1].n
    for p in patches:
        if p.task_type != 'rt_patch':
            gd.append(dict(
            left_edge = [a-b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
            right_edge = [a+b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
            level = p.level,
            dimensions=p.n,
            number_of_particles=0,
            density=p.data[:,:,:,p.varidx['d']],
            e_int = p.data[:,:,:,p.varidx['e']],
            dx=p.dx))

    return gd
    
def make_yt_dict_rt(patches):
    gd=[]
    n=patches[1].n
    for p in patches:
        if p.task_type == 'rt_patch':
            gd.append(dict(
            left_edge = [a-b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
            right_edge = [a+b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
            level = p.level,
            dimensions=p.n,
            number_of_particles=0,
            q = p.data[:,:,:,0]+1e-5))

    return gd    
    
def load_yt(dc,nx,sz=100,refine_by=3,geometry="cartesian",mass_unit=sc.scaling().m,length_unit=sc.scaling().l,
            time_unit=sc.scaling().t,velocity_unit=sc.scaling().v,periodicity=(False,False,False)):
    ds=yt.load_amr_grids(dc,nx,
                       bbox=np.array([[-sz,sz],[-sz,sz],[-sz,sz]]),
                       refine_by=3,
                       geometry="cartesian",
                       mass_unit=sc.scaling().m,
                       length_unit=sc.scaling().l,
                       velocity_unit=sc.scaling().v,
                       time_unit=sc.scaling().t,
                       periodicity=(False,False,False))
    return ds
