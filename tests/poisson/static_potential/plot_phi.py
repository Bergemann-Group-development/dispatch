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

    if x.ndim == 1 or y.ndim == 1:
        xx, yy = np.meshgrid(x,y,indexing='ij')
    else:
        xx, yy = x, y
    r = np.sqrt((xx-0.5)**2 + (yy-0.5)**2)
    phi_e = np.log((r / r0)**2 + 1.0)

    return phi_e

def rho_exact(x,y,r0=1.0):
    """
    \rho_e(r) = 4 * r_0^2 / (r^2 + r_0^2)^2

    """

    if x.ndim == 1 or y.ndim == 1:
        xx, yy = np.meshgrid(x,y,indexing='ij')
    else:
        xx, yy = x, y
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

    data = './data/input_1p_1280'
    iters = [0,10,50,90,150,190,250,290]
    pxsize = [-0.02,1.02]
    pysize = [-0.02,1.02]
    pltvar = ['d','rdd','phi','rdphi']
    cmin   = [ 1.0,  0.0,  0.0,   0.0]
    cmax   = [400., 1e-4,  3.9,   1.0]
    clog   = [   1,    0,    0,     0]
    showgrids = True
    showgridno = True
    ColMap = cm.RdBu_r
    r0 = 0.25

    L2_min, L2_max, L2_mean = [], [], []
    Linf_min, Linf_max, Linf_mean = [], [], []
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
            if j == 0:
                cax = cfig.add_subplot(gs[0,j])
            else:
                cax = cfig.add_subplot(gs[0,j],sharex=axs[0],sharey=axs[0])
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
                xp, yp, cut_mask = dg.plotting_grid(pp, x_i, 3, 0.0)
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
                    q0 = rho_exact(xp,yp,r0=r0)
                    var = relative_difference(q0[iindex,jindex,kindex].squeeze(),var.squeeze())
                elif pltvar[j] == 'rdphi':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['phi']]
                    q0 = phi_exact(xp,yp,r0=r0)
                    var = relative_difference(q0[iindex,jindex,kindex].squeeze(),var.squeeze())
                elif pltvar[j] == 'de':
                    var = rho_exact(xp,yp,r0=r0)
                    var = var[iindex,jindex]
                elif pltvar[j] == 'phie':
                    var = phi_exact(xp,yp,r0=r0)
                    var = var[iindex,jindex]
                else:
                    var = pp.data[iindex,jindex,kindex,pp.varidx[pltvar[j]]]

                print spltvar, ": ", var.min(), var.max(), var.shape
                xpp = xp[iindex,jindex,kindex].squeeze()
                ypp = yp[iindex,jindex,kindex].squeeze()

                if clog[j] == 1:
                    var = np.log10(var)
                    spltvar = r'log ' + spltvar

                conf = cax.contourf(xpp,ypp,var.squeeze(),levels,cmap=ColMap,
                             extend='both',zorder=-min(pp.dx))

                if (showgrids or pltvar[0] == 'gs'):
                    edges = dg.patch_grid_edges(pp,3,truecoords=True)
                    cax.add_patch(edges)
                # annotate patches with their ID (if desired)
                if showgridno:
                    offset = tuple([2*x for x in pp.nghost])
                    gn_xpos = xp[iindex,jindex,kindex][offset]
                    gn_ypos = yp[iindex,jindex,kindex][offset]
                    cax.annotate(str(pp.id),(gn_xpos,gn_ypos),color='k',weight='extra bold',
                                 fontsize=13,zorder=min(pp.dx))

            cbar = plt.colorbar(conf, ax=cax, fraction=0.03, pad=0.03)
            cbar.ax.annotate(spltvar,xy=(-0.10,-0.06),fontsize=14,xycoords='axes fraction')

        cfig.tight_layout()
        plt.draw()
        print "="*20

        # print L2 and Linf norm of the result
        L2_iter, Linf_iter = [], []
        for pp in patches:
            x_i = dg.dispatch_grid(pp)
            xp, yp, cut_mask = dg.plotting_grid(pp, x_i, 3, 0.0)
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

            phie = phi_exact(xp,yp,r0=r0)
            dphi = pp.data[:,:,:,pp.varidx['phi']] - phie
            L2d = np.linalg.norm(phie[iindex,jindex,kindex].flatten(),ord=2)
            L2n = np.linalg.norm(dphi[iindex,jindex,kindex].flatten(),ord=2)
            Linfd = np.linalg.norm(phie[iindex,jindex,kindex].flatten(),ord=np.inf)
            Linfn = np.linalg.norm(dphi[iindex,jindex,kindex].flatten(),ord=np.inf)
            L2_iter.append(L2n/L2d)
            Linf_iter.append(Linfn/Linfd)
            #print "i = {0}, p = {1}, n = {2}, L2-norm = {3:.3e}, Linf-norm = {4:.3e}".format(i,pp.id,pp.n,L2n/L2d,Linfn/Linfd)

        # accumulate L2, Linf results
        L2_min.append(np.min(L2_iter))
        L2_max.append(np.max(L2_iter))
        L2_mean.append(np.mean(L2_iter))
        Linf_min.append(np.min(Linf_iter))
        Linf_max.append(np.max(Linf_iter))
        Linf_mean.append(np.mean(Linf_iter))

    for i in xrange(len(iters)):
        print "i = {0}, iout = {1}".format(i,iters[i])
        print "L2 min = {0:.3e}, L2 max = {1:.3e}, L2 mean = {2:.3e}".format(L2_min[i],L2_max[i],L2_mean[i])
        print "Linf min = {0:.3e}, Linf max = {1:.3e}, Linf mean = {2:.3e}".format(Linf_min[i],Linf_max[i],Linf_mean[i])
        print "-"*40
