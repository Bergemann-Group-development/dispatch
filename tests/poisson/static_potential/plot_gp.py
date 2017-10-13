import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, patches as mplpatches, gridspec
import dispatch_utils as du
import dispatch_data as dd
import dispatch_grid as dg
from derived_variables import *
import itertools
import fnmatch



if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    data = './data/'
    iters = [2]
    pxsize = [ 0.0,1.0]
    pysize = [ 0.0,1.0]
    pltvar = ['d','phi','rd_d','rd_phi']
    cmin = [1.0, 0.1]
    cmax = [4.0, 0.4]
    clog = [  0,   0]
    showgrids = True
    ColMap = cm.RdBu_r

    for i in iters:
        patches = du.patches(data=data,iout=i)
        for pp in patches:
            pp.cache()
            pp.stripit()

        cfig = plt.figure("{0}".format(i),facecolor='w',figsize=(12.0,6.0))
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
