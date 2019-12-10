import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, patches as mplpatches, gridspec
import dispatch_utils as du
import dispatch_data as dd
import dispatch_grid as dg
from derived_variables import *
import itertools
import fnmatch
import os
from plot_phi import phi_exact, rho_exact

def amr_shading(levels=7):
    """
    Calculate what grey level is used to denote the level of refinement
    when plotting 1D AMR results.

    Plotting the level of refinement in this fashion is only "visually
    reasonable" if there are less than 10 levels, otherwise the contrast
    between the different grey levels becomes too small.

    Imported from soapp.plot.amr_common.

    >>> levels = 7
    >>> greylevs = amr_shading(levels)
    >>> print [g for g in greylevs]
    ['0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3']

    """

    if levels <= 9:
        dgrey = 0.1
    else:
        dgrey = 0.9 / levels
    
    # Note that a grey level of 1.0 (white) is reserved for the coarsest
    # level.
    fgreys = [1.0]
    for i in xrange(1,levels):
        fgreys.append(fgreys[i-1]-dgrey)

    return [str(i) for i in fgreys]

if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    data = './data/input_nested'
    iters = [0,1,3]
    pxsize = [ 0.0, 1.0]
    pltvar = ['phi','rdphi',   'd',  'rdd']
    cmin   = [  0.0, -0.020,   0.0,   -0.1]
    cmax   = [  5.0,  0.010,   4e2,    0.1]
    clog   = [    0,      0,     0,      0]
    slicedir = 1
    slicelocs = (0.460317,0.0)
    showgrids = False
    showgridno = True
    showghostzones = False
    skipped = []
    r0 = 0.1

    savecolr = dict()
    greys = amr_shading()
    scolrs = ['r','g','b','c','m','y','k','Brown','Coral','steelblue','limegreen','orange','fuchsia']
    scolrs = itertools.cycle(scolrs)
    istep_min, istep_max = [], []
    mi = 0
    for i in iters:
        if not os.path.exists(data):
            raise Exception("data dir not found!")
        patches = du.patches(data=data,iout=i)
        if len(patches) == 0: raise Exception("no patches found!")
        for pp in patches:
            pp.cache()
            pp.stripit()

        cfig = plt.figure("1-D {0}".format(i),facecolor='w',figsize=(5.0*len(pltvar),4.0))
        cfig.clf()
        gs = gridspec.GridSpec(1,len(pltvar))
        
        axs = []
        for j in xrange(len(pltvar)):
            if j == 0:
                cax = cfig.add_subplot(gs[0,j])
            else:
                cax = cfig.add_subplot(gs[0,j],sharex=axs[0])
            axs.append(cax)
            cax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            cax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            # ensure the axes remain visible by putting them at the "top"
            cax.set_zorder(99.0)

            spltvar = pltvar[j]
            if clog[j] == 1:
                spltvar = r'log ' + spltvar
                if cmin[j] != 0.0 and mi == 0: cmin[j] = np.log10(cmin[j])
                if cmax[j] != 0.0 and mi == 0: cmax[j] = np.log10(cmax[j])
            print cmin[j],cmax[j]
            cax.axis(xmin=pxsize[0],xmax=pxsize[1],ymin=cmin[j],ymax=cmax[j])

            for pp in patches:
                try:
                    skipped.index(pp.id)
                    print "Skipping ",pp.id," (user)"
                    continue
                except:
                    pass
                    
                if not savecolr.has_key(pp.id):
                    savecolr[pp.id] = scolrs.next()
                x_i = dg.dispatch_grid(pp)
                # ghost zones *can* now be included (or not)
                ione = (0,1)[pp.n[0] > 1]
                jone = (0,1)[pp.n[1] > 1]
                kone = (0,1)[pp.n[2] > 1]
                if np.any(pp.n != pp.nwrite):
                    iindex = slice(0,pp.n[0] + pp.nghost[0]*2)
                    jindex = slice(0,pp.n[1] + pp.nghost[1]*2)
                    kindex = slice(0,pp.n[2] + pp.nghost[2]*2)
                else:
                    iindex = slice(0,pp.n[0])
                    jindex = slice(0,pp.n[1])
                    kindex = slice(0,pp.n[2])
                if slicedir == 1:
                    if (slicelocs[0] < x_i[1].min() or slicelocs[0] > x_i[1].max()) and pp.n[1] > 1:
                        print "Skipping patch {0}...".format(pp.id)
                        continue
                    if (slicelocs[1] < x_i[2].min() or slicelocs[1] > x_i[2].max()) and pp.n[2] > 1:
                        print "Skipping patch {0}...".format(pp.id)
                        continue
                    dslice = np.abs(x_i[1] - slicelocs[0])
                    jindex, = np.where(dslice == dslice.min())[0]
                    dslice = np.abs(x_i[2] - slicelocs[1])
                    kindex, = np.where(dslice == dslice.min())[0]
                    if not showghostzones:
                        # don't show ghost zones in directions orthogonal to current slice
                        if jindex < pp.nghost[1] or jindex > (pp.n[1] + pp.nghost[1]):
                            print "Skipping patch {0}...".format(pp.id)
                            continue
                        if kindex < pp.nghost[2] or kindex > (pp.n[2] + pp.nghost[2]):
                            print "Skipping patch {0}...".format(pp.id)
                            continue
                    print "ID {0}: j index: {1}; k index: {2}".format(pp.id,jindex,kindex)
                    print x_i[1][jindex],x_i[2][kindex],jindex,kindex
                    cax.set_xlabel(r"$x$")

                xp = x_i[(slicedir-1)%2]
                yp = x_i[(slicedir-1+1)%2]
                zp = x_i[(slicedir-1+2)%2]

                print "Data:", pp.filename, "; step #: ",pp.istep
                if pltvar[j] == 'rdd':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['d']]
                    q0 = rho_exact(xp,yp,r0=r0)
                    var = (q0[iindex,jindex].squeeze() - var.squeeze()) / q0[iindex,jindex].squeeze()
                elif pltvar[j] == 'rdphi':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['phi']]
                    q0 = phi_exact(xp,yp,r0=r0)
                    var = (q0[iindex,jindex].squeeze() - var.squeeze()) / q0[iindex,jindex].squeeze()
                elif pltvar[j] == 'dphi':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['phi']]
                    q0 = phi_exact(xp,yp,r0=r0)
                    var = var.squeeze() - q0[iindex,jindex,kindex].squeeze()
                elif pltvar[j] == 'de':
                    var = rho_exact(xp,yp,r0=r0)
                    var = var[iindex,jindex]
                elif pltvar[j] == 'phie':
                    var = phi_exact(xp,yp,r0=r0)
                    var = var[iindex,jindex]
                elif pltvar[j] == 'dd':
                    var = pp.data[iindex,jindex,kindex,pp.varidx['d']]
                    q0 = rho_exact(xp,yp,r0=r0)
                    var = var.squeeze() - q0[iindex,jindex,kindex].squeeze()
                else:
                    var = pp.data[iindex,jindex,kindex,pp.varidx[pltvar[j]]]

                print spltvar, ": ", var.min(), var.max(), var.shape

                if clog[j] == 1:
                    var = np.log10(var)

                if pltvar[j] == 'phi':
                    q0 = phi_exact(xp,yp,r0=r0)
                    q0 = q0[iindex,jindex]
                    if clog[j] == 1:
                        q0 = np.log10(q0)
                    l0, = cax.plot(xp,q0.squeeze(),'k--',zorder=-min(pp.dx)*1.5)
                elif pltvar[j] == 'd':
                    q0 = rho_exact(xp,yp,r0=r0)
                    q0 = q0[iindex,jindex]
                    if clog[j] == 1:
                        q0 = np.log10(q0)
                    l0, = cax.plot(xp,q0.squeeze(),'k-',zorder=-min(pp.dx)*1.5)

                l, = cax.plot(xp,var.squeeze(),marker='o',color=savecolr[pp.id],ms=14.0/pp.level**1.25,
                              mew=0.0,zorder=-min(pp.dx),ls='none')
                cax.set_ylabel(spltvar)
                cax.grid(True)

                if (showgrids):
                    ylims = cax.get_ylim()
                    xregion = xp[pp.nghost[slicedir]:pp.n[slicedir]+pp.nghost[slicedir]]
                    rect = cax.fill_between(xregion,y1=ylims[0]-0.5*abs(ylims[0]),
                                            y2=ylims[1]+abs(ylims[1]),
                                            facecolor=greys[pp.level],zorder=-min(pp.dx)*2.5,
                                            edgecolor=greys[pp.level-1],lw=1.0,linestyle='dotted')
                    
                # annotate patches with their ID (if desired)
                if showgridno:
                    offset = 2*pp.nghost[slicedir]
                    gn_xpos = xp[offset]
                    cax.annotate(str(pp.id),(gn_xpos,0.8*cax.get_ylim()[1]),color='k',weight='extra bold',
                                 fontsize=13,zorder=min(pp.dx))

            print "-"*10

        cfig.tight_layout()
        cfig.subplots_adjust(bottom=0.14)
        plt.draw()
        print "="*20

        # increment counter (`iters` loop)
        mi += 1
