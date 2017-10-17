import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, patches as mplpatches, gridspec
import dispatch_utils as du
import dispatch_data as dd
import dispatch_grid as dg
from derived_variables import *
import itertools
import fnmatch

def phi_exact(x,y,r0=1.0):
    """
    \Phi_e(r) = ln( (r/r_0)^2 + 1)

    """

    xx, yy = np.meshgrid(x,y)
    r = np.sqrt((xx-0.5)**2 + (yy-0.5)**2)
    phi_e = np.log((r / r0)**2 + 1.0)

    return phi_e

def rho_exact(x,y,r0=1.0):
    """
    \rho_e(r) = 4 * r_0^2 / (r^2 + r_0^2)^2

    """

    xx, yy = np.meshgrid(x,y)
    r = np.sqrt((xx-0.5)**2 + (yy-0.5)**2)
    rho_e = 4.0 * r0**2 / (r**2 + r0**2)**2

    return rho_e

def relative_difference(qstd, qmeas):
    """Calculate the relative difference of a quantity."""
    
    return (np.abs(qstd - qmeas) / qstd)

if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    data = './data/'
    iters = [1,2]
    pxsize = [-0.02,1.02]
    pysize = [-0.02,1.02]
    pltvar = ['d','rdd','phi','rdphi']
    cmin   = [1.0,  0.0,  0.0,   0.0]
    cmax   = [16.,  1.0, 1.15,   1.0]
    showgrids = True
    ColMap = cm.RdBu_r
    r0 = 0.5

    for i in iters:
        patches = du.patches(data=data,iout=i)
        for pp in patches:
            pp.cache()
            pp.stripit()

        cfig = plt.figure("{0}".format(i),facecolor='w',figsize=(5.0*len(pltvar),4.0))
        cfig.clf()
        gs = gridspec.GridSpec(1,len(pltvar))
        
        axs = []
        for j in xrange(len(pltvar)):
            cax = plt.subplot(gs[0,j])
            axs.append(cax)
            cax.axis(xmin=pxsize[0],xmax=pxsize[1],ymin=pysize[0],ymax=pysize[1])
            cax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            cax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            # ensure the axes remain visible by putting them at the "top"
            cax.set_zorder(99.0)

            levels = np.linspace(cmin[j],cmax[j],num=50)
            spltvar = pltvar[j]

            for pp in patches:
                x_i = dg.dispatch_grid(pp)
                x, y, z = x_i
                # ghost zones *can* now be included (or not)
                ione = (0,1)[pp.n[0] > 1]
                jone = (0,1)[pp.n[1] > 1]
                kone = (0,1)[pp.n[2] > 1]
                if np.any(pp.n != pp.nwrite):
                    iindex = slice(pp.nghost[0] - ione,pp.n[0] + pp.nghost[0] + ione)
                    jindex = slice(pp.nghost[1] - jone,pp.n[1] + pp.nghost[1] + jone)
                    kindex = slice(pp.nghost[2] - kone,pp.n[2] + pp.nghost[2] + kone)
                else:
                    iindex = slice(0,pp.n[0])
                    jindex = slice(0,pp.n[1])
                    kindex = slice(0,pp.n[2])

                print "Data:", pp.filename, pp.istep
                if pltvar[j] == 'rdd':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['d']]
                    q0 = rho_exact(x,y,r0=r0)
                    print q0[iindex,jindex].shape,var.shape
                    var = relative_difference(q0[iindex,jindex],var.squeeze())
                elif pltvar[j] == 'rdphi':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['phi']]
                    q0 = phi_exact(x,y,r0=r0)
                    print q0[iindex,jindex].shape,var.shape
                    var = relative_difference(q0[iindex,jindex],var.squeeze())
                elif pltvar[j] == 'phie':
                    var = phi_exact(x,y,r0=r0)
                    var = var[iindex,jindex]
                else:
                    var = pp.data[iindex,jindex,kindex,pp.varidx[pltvar[j]]]

                print spltvar, ": ", var.min(), var.max(), var.shape
                conf = cax.contourf(x[iindex],y[jindex],var.squeeze(),levels,cmap=ColMap,
                             extend='both',zorder=-min(pp.dx))

                if (showgrids or pltvar[0] == 'gs'):
                    edges = dg.patch_grid_edges(pp,3,truecoords=True)
                    cax.add_patch(edges)
        
            cbar = plt.colorbar(conf, ax=cax, fraction=0.03, pad=0.03)
            cbar.ax.annotate(spltvar,xy=(-0.10,-0.06),fontsize=14,xycoords='axes fraction')

        cfig.tight_layout()
        plt.draw()
        print "="*20
