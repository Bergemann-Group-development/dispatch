from .atom import atom_class, CC
from .trans import trans_class
from .auxfuncs import vacuum2obs, continuum, freq_quad, convol
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import f90nml
from scipy.signal import argrelextrema
from scipy.io import FortranFile 
import os

cc = 2.99792458e10 # Speed of light [cm/s]

#=====================================================================
def prepare(self, dmax=None):
    #-----------------------------------------------------------------
    self.lam0 = vacuum2obs(np.array(self.alamb))
    self.lb = vacuum2obs(np.array(self.lam_blue))*(1-1e-7)
    self.lr = vacuum2obs(np.array(self.lam_red ))*(1+1e-7)
    self.i0 = np.searchsorted(self.run.lam, self.lam0)
    
    if dmax == None:
        q = freq_quad(self.nq, self.qmax, self.q0)
        lam = self.alamb / (q[:2]*self.run.qnorm/(CC*1e-5) + 1)
        self.dmax = np.abs(lam[1] - lam[0])
    else:
        self.dmax = dmax
    
    if self.run.nml['spectrum_params']['daa'] == 0:
        self.ib = self.i0 - np.argmax(self.run.dlam[self.i0::-1]>self.dmax) + 1
        self.ir = self.i0 + np.argmax(self.run.dlam[self.i0:]>self.dmax) + 1
        if self.ir == self.i0+1: self.ir = -1
    else:
        lb = self.run.nml['spectrum_params']['aa_blue']
        lr = self.run.nml['spectrum_params']['aa_red']
        if (self.lam0 < lr) and (self.lam0 > lb):
            self.ib = 0
            self.ir = self.run.lam.size
            
    self.lam = self.run.lam[self.ib:self.ir]
    self.flux = self.run.flux[self.ib:self.ir]
    self.ie = self.run.ie[:, self.ib:self.ir]
        
    self.prepared = True

#=====================================================================
def plot_trans(self, ang=None, imu=None, norm=True, dlam=None, \
               dmax=None, vdop=None, vshift=None, renorm=False, **kwargs):
    #-----------------------------------------------------------------
    if not self.prepared:
            self.prepare(dmax=dmax)
    
    if self.run.n_ang == 1: ang=0
        
    xx = self.lam
    
    if ang is not None:
        yy = self.ie[ang]
        plt.ylabel('Intensity µ={:.2f}'.format(self.run.vec[ang,-1]))
        
    elif imu is not None:
        yy, mu = self.get_Imu(imu, norm=False)
        plt.ylabel('Intensity µ={:.2f}'.format(mu))
        
    else:
        yy = self.flux
        plt.ylabel('Flux')
    
    
    if norm and not renorm:
        yy = yy / continuum(self.lam, yy)
        
    if vdop is not None or vshift is not None:
        sxx, syy = convol(xx, yy, vdop=vdop, vshift=vshift)
        yy = np.interp(xx, sxx, syy)

    if dlam==None:
        plt.xlim(self.run.lam[self.ib], self.run.lam[self.ir-1])
        plt.xlabel('Wavelength [Å]')
    elif dlam==True:
        xx = xx - self.lam0
        plt.xlabel('$\Delta\lambda$ [Å]')
    else:
        xx = xx - self.lam0
        plt.xlim(-dlam, dlam)
        plt.xlabel('$\Delta\lambda$ [Å]')
        
        if norm and renorm:
            mask = np.abs(xx) < dlam
            xx = xx[mask]
            yy = yy[mask]
            yy = yy / continuum(xx, yy)

    plt.title('Line = {:.2f} [Å], A({}) = {:.2f}'.format(self.lam0, self.el.id, self.el.abnd), fontsize='medium')
    plt.plot(xx, yy, **kwargs)
    plt.ticklabel_format(useOffset=False)

#=====================================================================
def get_Imu(self, imu=None, norm=True):
    #-----------------------------------------------------------------
    if not self.prepared: self.prepare()
    
    if imu is not None:
        mus, idx = np.unique(1-self.run.vec[:,-1], return_inverse=True)
        yy = np.mean(self.ie[idx==imu,:], axis=0)
        mu = 1-mus[imu]
        
    else:
        yy = self.flux
        mu = None
    
    if norm:
        yy = yy / continuum(self.lam, yy)
        
    return yy, mu
    
#=====================================================================
def calc_weq(self, ang=None):
    #-----------------------------------------------------------------
    if ang==None:
        yy = self.flux
    else:
        yy = self.ie[ang]
    
    yy = yy / continuum(self.lam, yy)
    weq = np.trapz(1 - yy, x=self.lam) * 1e3
    
    return weq

#=====================================================================
def calc_blueshift(self, imu=None):
    #-----------------------------------------------------------------
    nflux = self.get_Imu(imu=imu, norm=True)
    
    xx, indices = np.unique(self.lam.round(decimals=8), return_index=True)
    yy = nflux[indices]
    
    ex = argrelextrema(yy, np.less_equal)[0]
    iex = ex[np.argmin(np.abs(ex - self.i0 + self.ib))]
    imin3 = iex + np.array([-1, 0, 1])
    dep = 1 - yy[iex]
    
    p = np.polyfit(xx[imin3], yy[imin3], deg=2, rcond=1e-18)
    xmin = np.roots(np.polyder(p))[0]
    blue = (xmin / self.lam0 - 1) * cc * 1e-2
    
    return blue, dep


#=====================================================================
class read(object):
    #=================================================================
    def __init__(self, folder):
        #-------------------------------------------------------------
        folder = folder + os.sep
        self.nml = f90nml.read(folder + 'params.nml')
        
        df = pd.read_csv(folder + 'abund.txt', header=None, sep='\s+', na_filter=False)
        self.abundances = {k.lower().capitalize():v for k, v in df.values}
        
        if os.path.isfile(folder + 'atom.txt'):
            self.qnorm = self.nml['atom_params']['qnorm']
            self.atom = atom_class(folder + 'atom.txt', qnorm=self.qnorm)
            loggf_cut = self.nml['atom_params']['logf_cut']
            red_cut   = self.nml['atom_params']['red_cut']
            blue_cut  = self.nml['atom_params']['blue_cut']
            in_abnd   = self.nml['atom_params']['abundance']
            
            if self.atom.id in self.abundances:
                self.atom.abnd = self.abundances[self.atom.id]
            elif in_abnd > -99:
                self.atom.abnd = in_abnd
        
            if 'line_mask' in self.nml.keys():
                in_lines = self.nml['line_mask']['in_lines']
                lines = [self.atom.trans[i-1] for i in in_lines]
            else:
                lines = self.atom.trans[:self.atom.nline]
                
            lines = [line for line in lines if np.log10(line.fosc*line.lower.g) > loggf_cut]
            lines = [line for line in lines if line.lam_red > blue_cut and line.lam_blue < red_cut]
            
            self.nline = len(lines)
            self.line = lines
        
        print('Reading wavelengths')
        self.lam = np.fromfile(folder + 'lam.bin', dtype=np.float64)
        self.lam = vacuum2obs(self.lam[::-1])
        self.dlam = np.diff(self.lam)
        
        print('Reading Intensities')
        self.n_ang = len(glob.glob1(folder,"*_meta.txt"))
        self.vec = np.zeros([self.n_ang, 3])
        self.wts = np.zeros(self.n_ang)
        
        self.i3 = [None]*self.n_ang
        self.ie = [None]*self.n_ang
        
        self.nmu = np.unique(self.vec[:,-1]).size
        
        for i in range(self.n_ang):
            ray = 'ray{:03}'.format(i+1)
            meta = open(folder + ray + '_meta.txt')
            self.vec[i] = np.array(meta.readline().split()[1:], dtype=float)
            self.wts[i] = float(meta.readline().split()[1])
            info = meta.readline().split()
            file = info[0]
            prec = int(info[1])
            dim  = np.array(info[2:], dtype=int)
            
            if prec == 32: dtype=np.float32
            if prec == 64: dtype=np.float64
            
            self.i3[i] = np.memmap(folder + file, dtype=dtype).reshape(dim, order='F')[...,::-1]
            self.ie[i] = np.mean(self.i3[i], axis=(0,1))
            
        if (np.product(dim[:2]) > 1):
            self.dim = '3D'
        else:
            self.dim = '1D'
            
        self.i3 = np.array(self.i3)
        self.ie = np.array(self.ie)
        self.flux = np.sum(self.ie.T * self.wts, axis=-1) / np.sum(self.wts)
        self.flux_aa = self.flux*cc/(self.lam**2)*10
            
        if hasattr(self, 'atom'):
            #print('Prepairing transitions')
            for line in self.line:
                line.run = self
                line.prepared = False
            
        mesh = open(folder + 'atmos_mesh.txt', 'r')
        nx, ny, nz = [int(i) for i in mesh.readline().split()]
        dx, dy, dz = [float(i) for i in mesh.readline().split()]
        mesh.close()
        
        self.xx = np.linspace(0, (nx-1)*dx, nx)
        self.yy = np.linspace(0, (ny-1)*dy, ny)
        self.zz = np.linspace(0, (nz-1)*dz, nz)
        
        dat = FortranFile(folder + 'atmos_data.bin')
        self.temp = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        self.rho  = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        self.px   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        self.py   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        self.pz   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        self.pe   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        self.tau  = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
        dat.close()
        
        maxiter = self.nml['m3d_params']['maxiter']
        if (maxiter > 0):
            self.mode = 'NLTE'
        else:
            self.mode = 'LTE'
            
              
#=====================================================================
trans_class.prepare = prepare
trans_class.plot = plot_trans
trans_class.get_Imu = get_Imu
trans_class.calc_weq = calc_weq
trans_class.calc_blueshift = calc_blueshift