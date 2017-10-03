import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, patches as mplpatches
from dispatch_data import *
from dispatch_grid import *
from derived_variables import *
from patch_types import *
import itertools
import fnmatch

def dump_min_max(p):
    """
    Print the min and max of each variable in `p.data` (convenience function).

    """

    print "Time: {0}".format(p.time)
    for q,v in zip(p.varidx.iterkeys(),p.varidx.itervalues()):
        print "{0}: min = {1}; max = {2}".format(q,p.data[:,:,:,v].min(),p.data[:,:,:,v].max())

if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    gamma = 5.0/3.0
    datadir = './data/'
    prefix = ''
    ptimes = [0.0,0.10,0.20,0.35,0.45]
    t_eps = 1e-3
    pxsize = [-0.5,2.5]
    pltvar = [  'd', 'e1', 'v2', 'b2', 'et', 'v1']
    pvmin  = [  0.0,  0.0, -1.0, -0.1,  0.0, -0.1]
    pvmax  = [  1.1,  1.7,  0.1,  1.1,  2.6,  1.1]
    clog   = [    0,    0,    0,    0,    0,    0]
    #pltvar = [  'd', 'e1', 'v2', 'b2',  's', 'v1']
    #pvmin  = [  0.0,  0.0, -1.0, -0.1,  0.0, -0.1]
    #pvmax  = [  1.1,  1.7,  0.1,  1.1,  2.6,  1.1]
    #clog   = [    0,    0,    0,    0,    0,    0]
    plotdir = 1  # given with respect to the *root* grid (typically patch id = 1)
    slicelocs = (0.0,0.0) # coord. of slice to be taken in the other two directions
    nplots = 6
    nrows, ncols = 2,3
    showgrids = True
    output_png = False
    ThisPatch = ZeusMHDPatch
    ColMap = cm.RdBu_r
    verbose = False
    skipped = []

    xsize0 = 4.5
    ysize0 = 4.0

    m1 = 0 # counter 1
    root = None
    savecolr = dict()
    scolrs = ['r','g','b','c','m','y','k','Brown','Coral']
    scolrs = itertools.cycle(scolrs)

    if len(ptimes) > 30 and not output_png:
        print "Do you really want to open {0} windows?".format(len(ptimes))
        print "Try setting `output_png` to True."
        from sys import exit
        exit(1)

    dirs = os.listdir(datadir)
    dirs = [ff for ff in dirs if not ff.endswith('.log')]
    dirs = [ff for ff in dirs if not ff.endswith('.nq')]
    files = [dd+'/'+ff for dd in dirs for ff in fnmatch.filter(os.listdir(datadir+'/'+dd),'*.txt')]

    for pt in ptimes:

        cfig = plt.figure("t = {0:.3e}".format(pt),facecolor='w',figsize=(ncols*xsize0,nrows*ysize0))
        cfig.clf()
        if (nplots > nrows*ncols):
            raise Exception("Your asking for more plots than there are axes.")
        axs = []
        for i in range(1,nplots+1):
            if i == 1:
                cax = cfig.add_subplot(nrows,ncols,i)
            else:
                cax = cfig.add_subplot(nrows,ncols,i,sharex=axs[0])
            cax.axis(xmin=pxsize[0],xmax=pxsize[1])
            cax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            cax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            # ensure the axes remain visible by putting them at the "top"
            cax.set_zorder(99.0)
            axs.append(cax)

        patches = []
        times = []
        # in order to determine which patches have output at the desired time,
        # we need to look at ALL the patches! = P
        for f in files:
            ptemp = ThisPatch(datadir+f.strip('.txt'),verbose=verbose)
            if (abs(ptemp.time-pt)/(pt+1.0e-99) <= t_eps):
                if len(patches) > 0:
                    lastpt = patches[len(patches)-1]
                    if lastpt.id == ptemp.id:
                        if abs(lastpt.time-pt) > abs(ptemp.time-pt):
                            patches.pop(len(patches)-1)
                            patches.append(ptemp)
                        elif abs(lastpt.time-pt) <= abs(ptemp.time-pt):
                            # do nothing
                            pass
                    else:
                        patches.append(ptemp)
                else:
                    patches.append(ptemp)
            if (verbose): print '----'
        if len(patches) == 0:
            print "No patches found?! t = ",pt
            continue

        qmin = np.ones(nplots) * 1.0e99
        qmax = np.ones(nplots) *-1.0e99
        m2 = 0 # counter 2
        for p in patches:
            if p.kind == 'immersed_boundary_t':
                print "Skipping ",p.id," (immersed boundary)"
                continue
            try:
                skipped.index(p.id)
                print "Skipping ",p.id," (user)"
                continue
            except:
                pass
            if (root is None) or (root.id > p.id):
                root = p

            print "Data:", p.filename, p.time, pt
            times.append(p.time)

            x_i = dispatch_grid(p)
            gf = GeometricFactors()
            gf.init_grid(p, x_i)
            xp = x_i[plotdir-1]
            li = p.li[plotdir-1]-1
            ui = p.ui[plotdir-1]

            # ghost zones *can* now be included (or not)
            ione = (0,1)[p.n[0] > 1]
            jone = (0,1)[p.n[1] > 1]
            kone = (0,1)[p.n[2] > 1]
            if np.any(p.n != p.nwrite):
                iindex = slice(p.nghost[0] - ione,p.n[0] + p.nghost[0] + ione)
                jindex = slice(p.nghost[1] - jone,p.n[1] + p.nghost[1] + jone)
                kindex = slice(p.nghost[2] - kone,p.n[2] + p.nghost[2] + kone)
            else:
                iindex = slice(0,p.n[0])
                jindex = slice(0,p.n[1])
                kindex = slice(0,p.n[2])
            if plotdir == 1:
                xp = xp[iindex]

            # loop over axes/variables
            m3 = 0 # counter 3
            if not savecolr.has_key(p.id):
                savecolr[p.id] = scolrs.next()
            for ii in range(nplots):
                cax = axs[ii]
                cax.grid(True)

                if root.mesh_type == 'cylindrical':
                    if plotdir == 1:
                        cax.set_xlabel(r'$Z$')
                    elif plotdir == 2:
                        cax.set_xlabel(r'$R$')
                    elif plotdir == 3:
                        cax.set_xlabel(r'$\Phi$')
                elif root.mesh_type == 'spherical':
                    if plotdir == 1:
                        cax.set_xlabel(r'$r$')
                    elif plotdir == 2:
                        cax.set_xlabel(r'$\theta$')
                    elif plotdir == 3:
                        cax.set_xlabel(r'$\phi$')
                elif root.mesh_type == 'Cartesian':
                    if plotdir == 1:
                        cax.set_xlabel(r'$x$')
                    elif plotdir == 2:
                        cax.set_xlabel(r'$y$')
                    elif plotdir == 3:
                        cax.set_xlabel(r'$z$')

                if pltvar[ii] == 'v1':
                    var = p.data[iindex,jindex,kindex,p.varidx['p1']]
                    var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                elif pltvar[ii] == 'v2':
                    var = p.data[iindex,jindex,kindex,p.varidx['p2']]
                    var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                    var = var / gf['h2c'][iindex,np.newaxis,np.newaxis]
                    if p.mesh_type == 'spherical':
                        var = var / x_i[0][iindex,np.newaxis,np.newaxis]
                elif pltvar[ii] == 'v3':
                    var = p.data[iindex,jindex,kindex,p.varidx['p3']]
                    var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                    var = ( var / gf['h31c'][iindex,np.newaxis,np.newaxis] / gf['h32c'][np.newaxis,jindex,np.newaxis] )
                elif pltvar[ii] == 'divb':
                    var, sumd, normmax = divergence(p, (p.varidx['b1'],p.varidx['b2'],p.varidx['b3']), gf)
                    var = var[iindex,jindex,kindex]
                elif pltvar[ii] == 'e1':
                    if p.varidx.has_key('s'):
                        var = ( p.data[iindex,jindex,kindex,p.varidx['d']]**(gamma-1.0)
                                * np.exp(p.data[iindex,jindex,kindex,p.varidx['s']]*(gamma-1.0)
                                       / p.data[iindex,jindex,kindex,p.varidx['d']])
                                / (gamma-1.0) )
                    else:
                        var = p.data[iindex,jindex,kindex,p.varidx[pltvar[ii]]]
                else:
                    var = p.data[iindex,jindex,kindex,p.varidx[pltvar[ii]]]

                # record global max. and min. of data
                qmin[ii] = min(qmin[ii], np.nanmin(var))
                qmax[ii] = max(qmax[ii], np.nanmax(var))
    
                cpltvar = pltvar[ii]
                print cpltvar, ": ", var.min(), var.max(), var.shape
                if clog[ii] == 1:
                    var = np.log10(var)
                    cpltvar = r'log ' + cpltvar
                cax.set_ylabel(cpltvar)
                cax.axis(ymin=pvmin[ii],ymax=pvmax[ii])
                if m1 == 0 and m2 == 0 and clog[ii] == 1:
                    if cmin[ii] != 0.0: cmin[ii] = np.log10(cmin[ii])
                    if cmax[ii] != 0.0: cmax[ii] = np.log10(cmax[ii])
    
                l, = cax.plot(xp,var.squeeze(),marker='o',color=savecolr[p.id],zorder=-min(p.dx),mew=0.0,ls='None')
    
                if (showgrids or pltvar[0] == 'gs'):
                    cax.plot([x_i[plotdir-1][li],x_i[plotdir-1][li]],[-1.0e32,1.0e32],'k:',lw=2.0)
                    cax.plot([x_i[plotdir-1][ui],x_i[plotdir-1][ui]],[-1.0e32,1.0e32],'k:',lw=2.0)

                m3 += 1

            m2 += 1

        avg_time = np.mean(times)
        cfig.suptitle(r"$\langle t \rangle$ = {0:.3e}".format(avg_time),fontsize=16)
        cfig.tight_layout(rect=(0,0,1,0.97),pad=1.05)

        if output_png:
            currfile = "s{0:04d}.png".format(m1)
            print currfile
            cfig.savefig(currfile,dpi=100)
            plt.close(cfig)
        else:
            plt.draw()

        m1 += 1
        print "qmin, qmax: ", qmin, qmax
        print '================'
