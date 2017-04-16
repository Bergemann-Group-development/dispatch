import numpy as np
import matplotlib as mp
from dispatch_data import *
from dispatch_grid import *
from derived_variables import *
import itertools
import yt as yt
import scaling as sc

# a procedure to read in the data from patches. It uses dispatch_data.py to do that. Data from .txt files is
# stored in patches[] array, memory mapping is used for the actual data.
# input:
# snap - snapshot number
# datadir - self explanatory
# verbose - self explanatory
# ThisPatch - the type of patch we are dealing here with
# skip - patches to be skipped
# gamma - temorary here, later should be a memory entry
def get_data(snap,datadir,prefix='',verbose=False, ThisPatch='', skip=[], gamma = 1.4):
    ld = os.listdir(datadir)
    ld.sort()
    files = fnmatch.filter(ld,"*%.4i.txt"%snap)
    patches=[]
    for f in files:
        ptemp = ThisPatch(datadir+f.strip('.txt'),verbose=verbose)
        if (ptemp.id not in skip):
            patches.append(ptemp)
    print len(patches),"patches collected for snapshot", snap
    for p in patches:
        ul = np.zeros(3)
        ll=list(p.llc)
        for i in range(3):
            #ll[i] = p.pos[i] - p.dx[i]*0.5*p.n[i]
            ul[i] = p.llc[i] + p.dx[i]*float(p.n[i])
        #p.llc=tuple(ll)
        p.urc=tuple(ul)
        p.gamma = gamma                                                     
    return patches

def make_yt_dict(patches):
    gd=[]
    lbox=[0] * 3
    ubox=[0] * 3
    bbox=np.array([[0.,0.],[0.,0.],[0.,0.]])
    for p in patches:
        sg = (p.gamma-1.)*sc.scaling().temp
        if p.data.shape[1] > p.n[0]:
            ng=list(p.nghost)
            gn=list(p.gn)
            l=[0] *3
            u=[0] *3
            for cc in range(3):
                l[cc] = ng[cc]
                u[cc] = gn[cc] - ng[cc]
                if p.llc[cc] < lbox[cc]:
                    lbox[cc] = p.llc[cc]
                if p.urc[cc] > ubox[cc]:
                    ubox[cc] = p.urc[cc]
            gd.append(dict(
                left_edge = p.llc,
                right_edge = p.urc,
                level = p.level,
                dimensions=p.n,
                number_of_particles=0,
                density=p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']],
                e = p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['e']],
                e_d    = [a/b for a,b in zip(p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['e']],p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']])],
                dx=p.dx,
                velocity_x=[a/b for a,b in zip (p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['p1']],p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']])],
                velocity_y=[a/b for a,b in zip (p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['p2']],p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']])],
                velocity_z=[a/b for a,b in zip (p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['p3']],p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']])],
                temperature = [(a-0.5*((p1**2+p2**2+p3**2)/d))/d*sg for a,p1,p2,p3,d in zip(p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['e']], p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['p1']], p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['p2']], p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['p3']], p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],p.varidx['d']])],))
        else:
            for cc in range(3):
                if p.llc[cc] < lbox[cc]:
                    lbox[cc] = p.llc[cc]
                if p.urc[cc] > ubox[cc]:
                    ubox[cc] = p.urc[cc]
            gd.append(dict(
                left_edge = p.llc,
                right_edge = p.urc,
                level = p.level,
                dimensions=p.n,
                number_of_particles=0,
                density=p.data[:,:,:,p.varidx['d']],
                e_int = p.data[:,:,:,p.varidx['e']],
                e_int_d = [a/b for a,b in zip(p.data[:,:,:,p.varidx['e']],p.data[:,:,:,p.varidx['d']])],
                dx=p.dx,
                velocity_x=[a/b for a,b in zip (p.data[:,:,:,p.varidx['p1']],p.data[:,:,:,p.varidx['d']])],
                velocity_y=[a/b for a,b in zip (p.data[:,:,:,p.varidx['p2']],p.data[:,:,:,p.varidx['d']])],
                velocity_z=[a/b for a,b in zip (p.data[:,:,:,p.varidx['p3']],p.data[:,:,:,p.varidx['d']])],
                temperature = [(a-0.5*((p1**2+p2**2+p3**2)/d))/d*sg for a,p1,p2,p3,d in zip(p.data[:,:,:,p.varidx['e']], p.data[:,:,:,p.varidx['p1']], p.data[:,:,:,p.varidx['p2']], p.data[:,:,:,p.varidx['p3']], p.data[:,:,:,p.varidx['d']])]))
    for cc in range(3):
        bbox[cc][0]=lbox[cc]
        bbox[cc][1]=ubox[cc]
    patches[0].bbox=bbox
    return gd
    
    
    
def load_yt(dc,nx,bbox=np.array([[-100,100],[-100,100],[-100,100]]),refine_by=3,geometry="cartesian",time=0):
    ds=yt.load_amr_grids(dc,nx,
                       bbox=bbox,
                       refine_by=refine_by,
                       geometry="cartesian",
                       mass_unit=sc.scaling().m,
                       length_unit=sc.scaling().l,
                       velocity_unit=sc.scaling().v,
                       time_unit=sc.scaling().t,
                       sim_time = time,
                       periodicity=(False,False,False))
    return ds