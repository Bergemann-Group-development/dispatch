import numpy as np
import matplotlib as mp
from dispatch_data import *
from dispatch_grid import *
from derived_variables import *
import itertools
import yt as yt
import scaling as sc

def get_data(time,t_eps,datadir,prefix='',verbose=False, ThisPatch='', gamma=1.4):
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
    for p in patches:
        ul = np.zeros(3)
        for i in range(0,2):
            ll=list(p.llc)
            ll[i] = p.pos[i] - p.dx[i]*0.5*p.n[i]
            ul[i] = p.pos[i] + p.dx[i]*0.5*p.n[i]
            p.llc=tuple(ll)
            p.urc=tuple(ul)
        p.velx=np.zeros(shape=p.data.shape[0:3])
        p.vely=np.zeros(shape=p.data.shape[0:3])
        p.velz=np.zeros(shape=p.data.shape[0:3])
        p.temp=np.zeros(shape=p.data.shape[0:3])
        sg = (gamma-1.)*sc.scaling().temp
        for ii in range(0,p.data.shape[0]):
            for ij in range(0,p.data.shape[1]):
                for ik in range(0,p.data.shape[2]):
                    p.velx[ii,ij,ik] = p.data[ii,ij,ik,p.varidx['p1']]/p.data[ii,ij,ik,p.varidx['d']]
                    p.vely[ii,ij,ik] = p.data[ii,ij,ik,p.varidx['p2']]/p.data[ii,ij,ik,p.varidx['d']]
                    p.velz[ii,ij,ik] = p.data[ii,ij,ik,p.varidx['p3']]/p.data[ii,ij,ik,p.varidx['d']]
                    px = p.data[ii,ij,ik,p.varidx['p1']]
                    py = p.data[ii,ij,ik,p.varidx['p2']]
                    pz = p.data[ii,ij,ik,p.varidx['p3']]
                    k = 0.5*(px**2+py**2+pz**2)/p.data[ii,ij,ik,p.varidx['d']]
                    p.temp[ii,ij,ik] = (p.data[ii,ij,ik,p.varidx['e']] - k)/p.data[ii,ij,ik,p.varidx['d']] * sg                                                       
    return patches

def make_yt_dict(patches):
    gd=[]
    for p in patches:
        if p.data.shape[1] > p.n[0]:
            ng=list(p.nghost)
            gn=list(p.gn)
            l=ng
            u=[a-b for a,b in zip(gn,ng)]
            gd.append(dict(
                left_edge = [a-b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
                right_edge = [a+b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
                level = p.level-1,
                dimensions=p.n,
                number_of_particles=0,
                density=p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']],
                e = p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['e']],
                temperature= p.temp[l[0]:u[0],l[1]:u[1],l[2]:u[2]],
                e_d = [a/b for a,b in zip(p.data[:,:,:,p.varidx['e']],p.data[:,:,:,p.varidx['d']])],
                velocity_x=p.velx[l[0]:u[0],l[1]:u[1],l[2]:u[2]],
                velocity_y=p.vely[l[0]:u[0],l[1]:u[1],l[2]:u[2]],
                velocity_z=p.velz[l[0]:u[0],l[1]:u[1],l[2]:u[2]],
                dx=p.dx))
        else:
            gd.append(dict(
                left_edge = [a-b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
                right_edge = [a+b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
                level = p.level-1,
                dimensions=p.n,
                number_of_particles=0,
                density=p.data[:,:,:,p.varidx['d']],
                e_int = p.data[:,:,:,p.varidx['e']],
                e_int_d = [a/b for a,b in zip(p.data[:,:,:,p.varidx['e']],p.data[:,:,:,p.varidx['d']])],
                #velocity_x=p.velx,
                #velocity_y=p.vely,
                #velocity_z=p.velz,
                velocity_x=p.velx,
                velocity_y=p.vely,
                velocity_z=p.velz,
                dx=p.dx,
                temperature=p.temp))

    return gd
    
    
    
def load_yt(dc,nx,sz=100,refine_by=3,geometry="cartesian",mass_unit=sc.scaling().m,length_unit=sc.scaling().l,
            time_unit=sc.scaling().t,velocity_unit=sc.scaling().v,periodicity=(False,False,False),time=0):
    ds=yt.load_amr_grids(dc,nx,
                       bbox=np.array([[-sz,sz],[-sz,sz],[-sz,sz]]),
                       refine_by=3,
                       geometry="cartesian",
                       mass_unit=sc.scaling().m,
                       length_unit=sc.scaling().l,
                       velocity_unit=sc.scaling().v,
                       time_unit=sc.scaling().t,
                       sim_time = time,
                       periodicity=(False,False,False))
    return ds
