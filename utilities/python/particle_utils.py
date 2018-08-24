# -*- coding: utf-8 -*-
"""
    Particles object, with methods:

    part=particles(it=2,read=True)
    part.read_shell(r=100,dr=5)

Created on Wed Nov 29 09:25:50 2017

@author: Aake
"""
from __future__ import print_function
from scipy.io import FortranFile
import sys
if sys.version_info[0]!=2:
    from importlib import reload
from numpy import size,sqrt,product,mean,mod
import dispatch_utils; reload(dispatch_utils)
import dispatch_data; reload(dispatch_data)
from dispatch_utils import patch_files
from dispatch_data import Patch
from scaling import cgs, scaling
from time import time

class particles:
    def __init__(self,it=0,run='',data='../data',read=0,r=0,dr=1,gas_mass=0,verbose=0):
        ''' Initialize, by reading in the list of patch files '''
        start=time()
        self.files=patch_files(it,run,data)
        if gas_mass>0:
            self.gas_mass=gas_mass
        if verbose:
            print('{} files ({:.3f} sec)'.format(size(self.files),time()-start))
        if read:
            self.read_shell(r=r,dr=dr,verbose=verbose)

    def get_gas_mass(self,verbose=0):
        ''' Sum up the gas mass '''
        if verbose:
            start=time()
        print("computing gas mass ...")
        mass=0.0
        n=0
        for file in self.files:
            p=Patch(file)
            l=p.li
            u=p.ui
            mass+=product(p.size)*mean(p.data[l[0]:u[0],l[1]:u[1],l[2]:u[2],0])
            if verbose>1:
                n+=1
                if mod(n,200)==0: print(n,mass)
        self.gas_mass=mass
        if verbose:
            print('mass: {:.4f} ({:.2f} sec)'.format(mass,time()-start))

    def read_all(self,verbose=0):
        ''' Trigger reading of all patches (may take a long time) '''
        self.read_shell(r=0,verbose=verbose)

    def read_patch(self,file,verbose=0):
        ''' Read particles in a patch '''
        name=file+'.peb'
        f=FortranFile(name,'rb')
        fmt=f.readInts()
        if verbose>2:
            print('ioformat:',fmt[0])
        dim=f.readInts()
        dim=dim[0:2]
        a=f.readReals()
        a=a.reshape(dim[1],dim[0]).transpose()
        idv = f.readInts()
        idx = idv[0:dim[0]]
        f.close()
        dict={}
        for i in range(size(idx)):
            id=idx[i]
            dict[id]={'p':a[i,0:3],'v':a[i,0:3],'t':a[i,6],'s':a[i,7],'w':a[i,8]}
        return idx,dict

    def read_shell(self,r=100,dr=1,dust2gas=0.01,save=True,verbose=0):
        ''' Read particles in a shell, by default saving them '''
        dict={}
        n=0
        npatch=0
        rate=0.0
        start=time()
        t0=start
        if r==0:
            r2min=0.0
            r2max=1e30
        else:
            r2min=(r-dr)**2
            r2max=(r+dr)**2

        # open each patch file and check if relevant
        for file in self.files:
            p=Patch(file)
            rc=p.corner_radii()

            # If relevant, open the .peb file, and read the data
            if rc.min() < (r+dr) and rc.max() > (r-dr) or r==0:
                npatch+=1
                idx,dd=self.read_file(file)

                # If in the shell, sum up rate contribution, and add particle to dict
                if (r>0):
                    for i in range(size(idx)):
                        id=idx[i]
                        d=dd[id]
                        p=d['p']
                        p2=sum(p**2)
                        if p2>r2min and p2<r2max:
                            v=d['v']
                            w=d['w']
                            rate -= sum(p*v)/sqrt(p2)*w
                            n=n+1
                            if save:
                                dict[id]=d
                else:
                    for i in range(size(idx)):
                        id=idx[i]
                        dict[id]=dd[id]
                if verbose>1:
                    now=time()
                    print('{:.3f} sec'.format(now-start))
                    start=now
        if r==0:
            self.particles=dict
            print('{:.3f} sec'.format(time()-start))
            return
        if not hasattr(self,'gas_mass'):
            self.get_gas_mass(verbose=verbose)
        rate=rate/(2.0*dr)*self.gas_mass
        self.rate=rate
        if save:
            self.particles=dict
            if verbose>1:
                print('{} self.particles saved'.format(size(dict.keys())))
        if verbose:
            units=scaling()
            rate1=dust2gas*rate*cgs.yr/units.t*units.m/cgs.m_earth
            s='rate:{:9.2e} ={:9.2e} M_E/yr'.format(rate,rate1)
            s=s+', based on {} particles from {} patches'.format(n,npatch)
            print(s+' ({:.1f} sec)'.format(time()-t0))

    def read_plane(self,x=None,y=None,z=None,ds=1,wmin=1e-9,dust2gas=0.01,save=True,verbose=0):
        ''' Read particles in a plane, by default saving them '''
        dict={}
        start=time()

        def read_patch_plane(j,s,ds,c):
            ''' read particles in patch belonging to a plane '''
            if c[1] > (s-ds) and c[0] < (s+ds):
                if verbose>1:
                    print(file)
                idx,dd=self.read_patch(file)
                for i in range(size(idx)):
                    id=idx[i]
                    d=dd[id]
                    p=d['p']
                    if d['w']>wmin and p[j]>(s-ds) and p[j]<(s+ds):
                        dict[id]=d
                return 1
            else:
                return 0
        # open each patch file and check if relevant
        npatch=0
        for file in self.files:
            p=Patch(file)
            c=p.active_corners
            if type(x) != type(None):
                npatch+=read_patch_plane(0,x,ds,c[:,0])
            elif type(y) != type(None):
                npatch+=read_patch_plane(1,y,ds,c[:,1])
            elif type(z) != type(None):
                npatch+=read_patch_plane(2,z,ds,c[:,2])
            if verbose>1:
                now=time()
                print('{}: {:.3f} sec'.format(file,now-start))
                start=now
        self.particles=dict
        if verbose>0:
                print('{} particles from {} patches ({:.2f} sec)'.format(size(dict.keys()),npatch,time()-start))

    def shell(self,r=100,dr=1,dust2gas=0.01,verbose=0):
        ''' evaluate rate of particles in a shell '''
        n=0
        rate=0.0
        t0=time()
        if r==0:
            r2min=0.0
            r2max=1e30
        else:
            r2min=(r-dr)**2
            r2max=(r+dr)**2
        for id in (self.particles).keys():
            part=self.particles[id]
            p=part['p']
            p2=sum(p**2)
            if p2>r2min and p2<r2max:
                v=part['v']
                w=part['w']
                rate -= sum(p*v)/sqrt(p2)*w
                n=n+1
        if not hasattr(self,'gas_mass'):
            self.get_gas_mass(verbose=verbose)
        rate=rate/(2.*dr)*self.gas_mass
        self.rate=rate
        if verbose:
            units=scaling()
            rate1=dust2gas*rate*cgs.yr/units.t*units.m/cgs.m_earth
            s='rate:{:9.2e} ={:9.2e} M_E/yr, based on {} particles.'.format(rate,rate1,n)
            s=s+' Time used: {:.1f} sec.'.format(time()-t0)
            print(s)
