import numpy as np
import matplotlib as mp
from dispatch_data import *
from dispatch_grid import *
from derived_variables import *
import itertools
import yt as yt
import scaling as sc

ptimes=[0.2]
t_eps=0.01
ThisPatch = RamsesMHDPatch
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
        ul=np.zeros(3)
        for i in range(0,2):
            ll=list(p.llc)
            ll[i] = p.pos[i] - p.dx[i]*0.5*p.n[i]
            ul[i] = p.pos[i] + p.dx[i]*0.5*p.n[i]
            p.llc=tuple(ll)
            p.urc=tuple(ul)
        #p.velx0 = np.zeros([p.n[0],p.n[1],p.n[2]])
        #p.vely0 = np.zeros([p.n[0],p.n[1],p.n[2]])
        #p.velz0 = np.zeros([p.n[0],p.n[1],p.n[2]])
        #for ii in range(0,p.n[0]):
         #   for ij in range(0,p.n[1]):
          #      y=p.llc[1]+ij*p.dx[1]-pl_pos[1]
           #     for ik in range(0,p.n[2]):
            #        z=p.llc[2]+ik*p.dx[2]-pl_pos[2]
                    #'''
                    # A rotation profile, where the rotation center is at y=+1 AU, and
                    # the direction of rotation is +x, so the head wind is in the -x
                    # direction.  The slower rotation at larger radii means that Ux is
                    # proportional to +y.  
                    #'''
                    #dp1=rho_gas*(0.5*omega[2]*y) * rho_gas
                    #dp2=v_accretion*rho_gas
                    #dp3=0.0

                    #'''
                    #the vertical and radial rotating coordinate frame accelerations.
                    #The force in the radial direction increases with increasing y,
                    ##giving an equilbrium velocity in the x-direction that is positive
                    #(faster rotation) for a positive y-displacement.  Vertically, there
                    #is a net force that comes from the projection of the force of
                    #gravity onto the z-axis, and has the same scaling constant as the
                    #radial force, but the opposite sign.
                    #'''
                    #dp2=dp2 + 3.0*omega[2]**2*y
                    #dp3=dp3 -     omega[2]**2*z

                    #'''
                    #the classical Coriolis force term: -2*(\Omega x U). A positive omega(3)
                    #is motion in the positive x-direction, if y points towards the system center.
                    #Hence, a positive velocity in the x-direction increases the centripetal force,
                    #and leads to acceleration in the negative y-direction -- vice versa for
                    #motion in the y-direction; a motion inwards (+y) leads to acceleration in
                    #the forward (+x) direction
                    #'''
                    #ddp1=dp1*rho_gas-2.*(omega[1]*dp3-omega[2]*dp2)
                    #ddp2=dp2*rho_gas-2.*(omega[2]*dp1-omega[0]*dp3)
                    #ddp3=dp3*rho_gas-2.*(omega[0]*dp2-omega[1]*dp1)

                    #p.velx0[ii,ij,ik]=ddp1/rho_gas
                    #p.vely0[ii,ij,ik]=ddp2/rho_gas
                    #p.velz0[ii,ij,ik]=ddp3/rho_gas


        p.velx=get_velocity(p.n,p.data[:,:,:,p.varidx['d']],p.data[:,:,:,p.varidx['p1']])
        p.vely=get_velocity(p.n,p.data[:,:,:,p.varidx['d']],p.data[:,:,:,p.varidx['p2']])
        p.velz=get_velocity(p.n,p.data[:,:,:,p.varidx['d']],p.data[:,:,:,p.varidx['p3']])
        
    return patches

def make_yt_dict(patches):
    gd=[]
    n=patches[1].n
    for p in patches:
        gd.append(dict(
            left_edge = [a-b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
            right_edge = [a+b for a,b in zip(p.pos,[c*0.5*d for c,d in zip(p.dx,p.n)])],
            level = p.level,
            dimensions=p.n,
            number_of_particles=0,
            density=p.data[:,:,:,p.varidx['d']],
            e_int = p.data[:,:,:,p.varidx['e']],
            e_int_d = [a/b for a,b in zip(p.data[:,:,:,p.varidx['e']],p.data[:,:,:,p.varidx['d']])],
            velocity_x=p.velx,#[a-b for a,b in zip(p.velx,p.velx0)], #p.velx
            velocity_y=p.vely,#[a-b for a,b in zip(p.velx,p.velx0)], #p.vely,
            velocity_z=p.velz,#[a-b for a,b in zip(p.velx,p.velx0)], #p.velz,
         #   velx0=p.velx0,
         #   vely0=p.vely0,
         #   velz0=p.velz0,
            dx=p.dx))

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
