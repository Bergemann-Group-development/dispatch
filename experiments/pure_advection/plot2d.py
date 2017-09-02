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

    gamma = 1.4
    datadir = './data/'
    prefix = ''
    ptimes = [0.0, 0.5, 1.0, 1.5, 2.0]
    #ptimes = 0.629*np.arange(0,2)
    t_eps = 1e-4
    pxsize = [-0.51,0.51]
    pysize = [-0.51,0.51]
    pltvar = ['d']
    cmin   = [ 0.0]
    cmax   = [ 4.2]
    clog   = [   0]
    slicedir = 3  # given with respect to the *root* grid (typically patch id = 1)
    sliceloc = 0.0 # coord. of slice to be taken
    showgrids = True
    showgridno = True
    output_png = False
    truecoords = False
    ThisPatch = ZeusMHDPatch
    ColMap = cm.RdBu_r
    verbose = False
    skipped = []

    m1 = 0 # counter 1
    root = None
    if len(ptimes) > 30 and not output_png:
        print "Do you really want to open {0} windows?".format(len(ptimes))
        print "Try setting `output_png` to True."
        from sys import exit
        exit(1)
    for pt in ptimes:

        cfig = plt.figure("t = {0:.3e}".format(pt),facecolor='w',figsize=(7.0,6.0))
        cfig.clf()
        cax = cfig.add_subplot(1,1,1)
        cax.axis(xmin=pxsize[0],xmax=pxsize[1],ymin=pysize[0],ymax=pysize[1])
        #cax.axes.set_aspect('equal')
        cax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        cax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        # ensure the axes remain visible by putting them at the "top"
        cax.set_zorder(99.0)

        files = fnmatch.filter(os.listdir(datadir), prefix + "*.txt")
        patches = []
        times = []
        # in order to determine which patches have output at the desired time,
        # we need to read ALL the patches! = P
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

        qmin, qmax = 1.0e99, -1.0e99
        m2 = 0 # counter 2
        for p in patches:
            if p.kind == 'immersed_boundary_t':
                circle = mplpatches.Circle((p.pos[0],p.pos[1]),radius=pradius,color='k')
                cax.add_patch(circle)
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
            if truecoords:
                if root.mesh_type == 'cylindrical':
                    if slicedir == 1:
                        cax.set_xlabel(r'$r$')
                        cax.set_ylabel(r'$\phi$')
                    elif slicedir == 2:
                        cax.set_xlabel(r'$z$')
                        cax.set_ylabel(r'$\phi$')
                    elif slicedir == 3:
                        cax.set_xlabel(r'$r$')
                        cax.set_ylabel(r'$z$')
                elif root.mesh_type == 'spherical':
                    if slicedir == 1:
                        cax.set_xlabel(r'$\theta$')
                        cax.set_ylabel(r'$\phi$')
                    elif slicedir == 2:
                        cax.set_xlabel(r'$r$')
                        cax.set_ylabel(r'$\phi$')
                    elif slicedir == 3:
                        cax.set_xlabel(r'$r$')
                        cax.set_ylabel(r'$\theta$')
            else:
                if slicedir == 1:
                    cax.set_xlabel(r'$y$')
                    cax.set_ylabel(r'$z$')
                elif slicedir == 2:
                    cax.set_xlabel(r'$x$')
                    cax.set_ylabel(r'$z$')
                elif slicedir == 3:
                    cax.set_xlabel(r'$x$')
                    cax.set_ylabel(r'$y$')

            print "Data:", p.filename, p.time, pt, p.pos[0], p.llc[0]
            times.append(p.time)

            x_i = dispatch_grid(p,truecoords=truecoords)
            gf = GeometricFactors()
            gf.init_grid(p, x_i)
            xp, yp, cut_mask = plotting_grid(root, p, x_i, slicedir, sliceloc,
                                             truecoords=truecoords,sliceeps=0.1)
            print x_i[0][0]
            # ghost zones *can* now be included (or not)
            ione = (0,1)[p.n[0] > 1]
            jone = (0,1)[p.n[1] > 1]
            kone = (0,1)[p.n[2] > 1]
            if np.any(p.n != p.nwrite):
                iindex = slice(p.nghost[0],p.n[0] + p.nghost[0])
                jindex = slice(p.nghost[1],p.n[1] + p.nghost[1])
                kindex = slice(p.nghost[2],p.n[2] + p.nghost[2])
            else:
                iindex = slice(0,p.n[0])
                jindex = slice(0,p.n[1])
                kindex = slice(0,p.n[2])

            if pltvar[0] == 'v1':
                var = p.data[iindex,jindex,kindex,p.varidx['s1']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
            elif pltvar[0] == 'v2':
                var = p.data[iindex,jindex,kindex,p.varidx['s2']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                var = var / gf['h2c'][iindex,np.newaxis,np.newaxis]
                if p.mesh_type == 'spherical':
                    var = var / x_i[0][iindex,np.newaxis,np.newaxis]
            elif pltvar[0] == 'v3':
                var = p.data[iindex,jindex,kindex,p.varidx['s3']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                var = ( var / gf['h31c'][iindex,np.newaxis,np.newaxis] / gf['h32c'][np.newaxis,jindex,np.newaxis] )
            elif pltvar[0] == 'divb':
                var, sumd, normmax = divergence(p, (p.varidx['b1'],p.varidx['b2'],p.varidx['b3']), gf)
                var = var[iindex,jindex,kindex]
            else:
                var = p.data[iindex,jindex,kindex,p.varidx[pltvar[0]]]

            # record global max. and min. of data
            qmin = min(qmin, np.nanmin(var))
            qmax = max(qmax, np.nanmax(var))

            # mask data for cutting
            var = np.ma.array(var,mask=cut_mask)
            cpltvar = pltvar[0]
            print cpltvar, ": ", var.min(), var.max(), var.shape
            # Note that these values, in general, will differ between patches.
            if cmin[0] == 0.0 and cmax[0] == 0.0:
                cmin[0] = qmin
                cmax[0] = qmax
            if clog[0] == 1:
                var = np.log10(var)
                cpltvar = r'log ' + cpltvar
            if m1 == 0 and m2 == 0 and clog[0] == 1:
                if cmin[0] != 0.0: cmin[0] = np.log10(cmin[0])
                if cmax[0] != 0.0: cmax[0] = np.log10(cmax[0])

            levels = np.linspace(cmin[0], cmax[0], num=50)
            conf = cax.contourf(xp[iindex,jindex,kindex].squeeze(),yp[iindex,jindex,kindex].squeeze(),
                                var.squeeze(),levels,cmap=ColMap,extend='both',zorder=p.quality*0.1)

            if (showgrids or pltvar[0] == 'gs'):
                edges = patch_grid_edges(p,slicedir,truecoords=truecoords)
                cax.add_patch(edges)
                #for c in conf.collections:
                #    c.set_clip_path(edges)

            # annotate patches with their ID (if desired)
            if p.id != 1 and showgridno:
                offset = tuple([3*x for x in p.nghost])
                gn_xpos = xp[iindex,jindex,kindex][offset]
                gn_ypos = yp[iindex,jindex,kindex][offset]
                cax.annotate(str(p.id),(gn_xpos,gn_ypos),color='k',weight='extra bold',
                            fontsize=13*(root.level+1)/(p.level+1),zorder=min(p.dx)*10.0)

            m2 += 1

        cbar = plt.colorbar(conf, ax=cax, fraction=0.03, pad=0.03)
        cbar.ax.tick_params(axis='both',which='major',labelsize=16)
        cbar.ax.annotate(cpltvar,xy=(-0.10,-0.06),fontsize=14,xycoords='axes fraction')

        avg_time = np.mean(times)
        cfig.suptitle(r"$\langle t \rangle$ = {0:.3e}".format(avg_time),fontsize=16)
        cfig.tight_layout()

        if output_png:
            currfile = "p{0:04d}.png".format(m1)
            print currfile
            cfig.savefig(currfile,dpi=100)
            plt.close(cfig)
        else:
            plt.draw()

        m1 += 1
        print "qmin, qmax: ", qmin, qmax
        print '================'
