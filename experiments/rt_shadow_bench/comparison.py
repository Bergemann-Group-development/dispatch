# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import ticker, cm, pyplot as plt
import constants as cnst

"""
Provide a simple estimation of the heating rate along the x-axis (assuming
y = z = 0) for the RT shadow benchmark.

"""

def density (x, d0=1.0, d1=10.0, a=0.10):
    """Calculate the density as a function of x."""

    r = (x / a)**2
    d = d0 + (d1 - d0) / (1 + np.exp(10.0 * (r - 1.0)))

    return d

def extinction (d, T, alpha0 = 1.0, T0=290.0, d0=1.0):
    """Calculate the opacity, \kappa, given the temperature and density."""

    alpha = alpha0 * (T / T0)**(-3.5) * (d / d0)**2
    alpha /= d

    return alpha

def source (T):
    """Calculate the blackbody source function for a given temperature."""

    return cnst.SBc * T**4

if __name__ == "__main__":
    # zc = zone-centred; fc = face-centred

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    make_plots = True
    ncells = 256

    ds = 1./ncells
    d0 = 1.0
    T0 = 290.0
    d1 = 10.0 * d0
    T1 = d0 * T0 / d1
    Thot = 6.0 * T0

    x = np.linspace(-0.5, 0.5, ncells+1)
    d = density(x, d0=d0, d1=d1)
    T = d0 * T0 / d
    alpha = extinction(d, T) # zc
    dtau = alpha * ds # zc

    # Approach #1: Divergence of radiative flux
    Flux_in = source(Thot)
    # note that the flux is staggered back by a half-zone
    Flux = np.empty(d.size+1)
    Flux[0] = Flux_in
    Flux[1:] = Flux_in * np.exp(-np.cumsum(dtau)) # fc
    Q1 = -np.diff(Flux) / ds
    
    # Approach #2: Integral solution of Q+ = I+ - S
    Src_in = source(Thot) - source(T0)
    Src = source(T) # zc
    ai = np.exp(-dtau) # zc
    bi = 1.0 - ai # zc
    ci = bi - dtau * ai #zc
    dtau_f = np.empty_like(Src)
    for i in xrange(1,ncells):
        dtau_f[i] = 0.5 * (dtau[i] + dtau[i-1])
    dtau_c = np.empty_like(Src)
    for i in xrange(0,ncells-1):
        dtau_c[i] = 0.5 * (dtau_f[i] + dtau_f[i+1])
    dSdtau = np.diff(Src) / dtau_f[:ncells] #fc
    Sprime = np.empty_like(Src)
    Sprimeprime = np.empty_like(Src)
    Qplus = np.empty_like(Src)
    for i in xrange(1,ncells-1):
        Sprime[i] = 0.5 * (dSdtau[i] * dtau_f[i+1] + dSdtau[i+1] * dtau_f[i] ) / dtau_c[i] # zc
        Sprimeprime[i] = (dSdtau[i+1] - dSdtau[i]) / dtau_c[i] # zc
    Qplus[0] = Src_in
    for i in xrange(1,ncells):
        Qplus[i] = ai[i] * Qplus[i-1] - bi[i] * Sprime[i] + ci[i] * Sprimeprime[i] # zc
    # the actual heating rate, with units energy per unit volume per unit time
    Q2 = alpha * Qplus

    # Relative difference between heating rates calculating with the two methods
    rd = np.abs(Q1 - Q2) / np.abs(Q2)
    print "Max Q1: {0:.5e}".format(np.nanmax(Q1[:ncells]))
    print "Max Q2: {0:.5e}".format(np.nanmax(Q2[:ncells]))
    print "Max relative difference: {0:.5e}".format(np.nanmax(rd))

    # Make some plots.
    if make_plots:
        fig = plt.figure('heating rates',facecolor='w',figsize=(14.0,8.0))
        fig.clf()
        ax1 = fig.add_subplot(2,3,1)
        ax2 = fig.add_subplot(2,3,2)
        ax3 = fig.add_subplot(2,3,3)
        ax4 = fig.add_subplot(2,3,4)
        ax5 = fig.add_subplot(2,3,5)

        ax1.semilogy(x,Src,'ko',mew=0.0)
        ax1.set_ylabel(r"$S = \sigma T^4$")
        ax2.semilogy(x,dtau,'yo',mew=0.0)
        ax2.set_ylabel(r"$\Delta \tau$")
        ax3.plot(x[:ncells],dSdtau,'co',mew=0.0)
        ax3.set_ylabel(r"$dS/d\tau$")
        ax4.plot(x,Q1,'bo',mew=0.0,label=r"$\nabla\cdot\vec{F}$")
        ax4.plot(x,Q2,'go',mew=0.0,label=r"Integral")
        ax4.set_ylabel(r"$Q$")
        ax4.legend(numpoints=1,fontsize=12,frameon=False,handletextpad=0.0)
        ax4.set_xlabel(r"$x$")
        ax5.semilogy(x,rd,'ro',mew=0.0)
        ax5.set_ylabel("Rel. Diff.")
        ax5.set_xlabel(r"$x$")
        ax5.axis(ymax=3.0e3)

        fig.tight_layout()
        plt.show()
