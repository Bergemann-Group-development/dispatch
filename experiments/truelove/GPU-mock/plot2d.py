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

def set_variable_name(v, passthrough=False):
    """
    Given a variable name, convert it to its "pretty-print" form (e.g., Greek
    letter) if `passthrough` is `False`.
    Variables must be added by hand.

    """

    ppv = v
    if v == 'd':
        ppv = r"$\rho$"

    return ppv

if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    gamma = 1.0001
    datadir = './data/'
    prefix = ''
    ptimes = [0.635]
    t_eps = 1e-3
    #pxsize = [-1.0,1.0]
    #pysize = [-1.0,1.0]
    pxsize = [-0.029,0.029]
    pysize = [-0.029,0.029]
    #xscaling = 3.12e16
    #xscaling = 2.89e17
    xscaling = 2.0e17
    #xscaling = None
    xunits = r"cm"
    yunits = r"cm"
    pltvar = ['d']
    cmin   = [2e3]
    cmax   = [1.75e7]
    clog   = [  1]
    vecvar = ['v1','v2']
    #vecvar = None
    ivfreq = 1./81. # at the level of the root patch
    refine_ratio = 3
    #vscaling = [6.549e-17]
    #vscaling = [8.240e-20]
    vscaling = [2.4861e-19]
    #vscaling = [None]
    vunits = [r"g cm$^{-3}$"]
    #vunits = None
    #vecscaling = 7.123e7
    #vecscaling = 7.598e4
    vecscaling = 9.1324e4
    #vecscaling = 1.0
    vecunits = r"cm s$^{-1}$"
    vpscale = 30.0
    global_llc = [-1.0,-1.0,-1.0]
    slicedir = 3  # given with respect to the *root* grid (typically patch id = 1)
    sliceloc = 0.0 # coord. of slice to be taken
    sliceeps = 1.0e-4
    showgrids = True
    showgridno = False
    output_png = False
    truecoords = True
    ThisPatch = ZeusMHDPatch
    ColMap = cm.RdBu_r
    verbose = False
    suppress_title = True
    skipped = [1,32,69,93,96,99,120,339,342,282,285,288,309]
    #skipped = [1,19,46,73,69,64,82,76,70,43,55,49,37] # skip patches entirely covered by others (manual for now)

    m1 = 0 # counter 1
    root = None
    if len(ptimes) > 30 and not output_png:
        print "Do you really want to open {0} windows?".format(len(ptimes))
        print "Try setting `output_png` to True."
        from sys import exit
        exit(1)

    dirs = os.listdir(datadir)
    dirs = [ff for ff in dirs if not ff.endswith('.log')]
    dirs = [ff for ff in dirs if not ff.endswith('.nq')]
    files = [dd+'/'+ff for dd in dirs for ff in fnmatch.filter(os.listdir(datadir+'/'+dd),'*.txt')]

    if xscaling is not None:
        pxsize = [px * xscaling for px in pxsize]
        pysize = [py * xscaling for py in pysize]
    sxunits, syunits = "", ""
    if xscaling is not None and xunits is not None:
        sxunits = " [" + xunits + "]"
    if xscaling is not None and yunits is not None:
        syunits = " [" + yunits + "]"

    for pt in ptimes:

        cfig = plt.figure("t = {0:.3e}".format(pt),facecolor='w',figsize=(8.5,7.5))
        cfig.clf()
        cax = cfig.add_subplot(1,1,1)
        cax.axis(xmin=pxsize[0],xmax=pxsize[1],ymin=pysize[0],ymax=pysize[1])
        #cax.axes.set_aspect('equal')
        cax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        cax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        # ensure the axes remain visible by putting them at the "top"
        cax.set_zorder(99.0)

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
        vmagmax_global = 0.0
        m2 = 0 # counter 2
        pused = []
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
            pused.append(p)
            if truecoords:
                if root.mesh_type == 'cylindrical':
                    if slicedir == 1:
                        cax.set_xlabel(r'$r$' + sxunits)
                        cax.set_ylabel(r'$\phi$' + syunits)
                    elif slicedir == 2:
                        cax.set_xlabel(r'$z$' + sxunits)
                        cax.set_ylabel(r'$\phi$' + syunits)
                    elif slicedir == 3:
                        cax.set_xlabel(r'$r$' + sxunits)
                        cax.set_ylabel(r'$z$' + syunits)
                elif root.mesh_type == 'spherical':
                    if slicedir == 1:
                        cax.set_xlabel(r'$\theta$' + sxunits)
                        cax.set_ylabel(r'$\phi$' + syunits)
                    elif slicedir == 2:
                        cax.set_xlabel(r'$r$' + sxunits)
                        cax.set_ylabel(r'$\phi$' + syunits)
                    elif slicedir == 3:
                        cax.set_xlabel(r'$r$' + sxunits)
                        cax.set_ylabel(r'$\theta$' + syunits)
                elif root.mesh_type == 'Cartesian':
                    if slicedir == 1:
                        cax.set_xlabel(r'$y$' + sxunits)
                        cax.set_ylabel(r'$z$' + syunits)
                    elif slicedir == 2:
                        cax.set_xlabel(r'$x$' + sxunits)
                        cax.set_ylabel(r'$z$' + syunits)
                    elif slicedir == 3:
                        cax.set_xlabel(r'$x$' + sxunits)
                        cax.set_ylabel(r'$y$' + syunits)
            else:
                if slicedir == 1:
                    cax.set_xlabel(r'$y$' + sxunits)
                    cax.set_ylabel(r'$z$' + syunits)
                elif slicedir == 2:
                    cax.set_xlabel(r'$x$' + sxunits)
                    cax.set_ylabel(r'$z$' + syunits)
                elif slicedir == 3:
                    cax.set_xlabel(r'$x$' + sxunits)
                    cax.set_ylabel(r'$y$' + syunits)

            print "Data:", p.filename, p.time, pt

            x_i = dispatch_grid(p)
            gf = GeometricFactors()
            gf.init_grid(p, x_i)
            xp, yp, cut_mask = plotting_grid(p, x_i, slicedir, sliceloc,
                                             truecoords=truecoords,sliceeps=sliceeps)
            ivfreql = max(int(ivfreq * refine_ratio**p.level), 1)
            if xscaling is not None:
                xp *= xscaling
                yp *= xscaling
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

            act_slice = 0.0
            if slicedir == 1 and ione > 0:
                diff = np.abs(x_i[0] - sliceloc)
                iindex, = np.where(diff < sliceeps)
                iindex = iindex[np.logical_and(iindex >= p.li[0], iindex < p.ui[0])]
                act_slice = x_i[0][iindex.min()]
                iindex  = slice(iindex.min(),iindex.max()+1)
            elif slicedir == 2 and jone > 0:
                diff = np.abs(x_i[1] - sliceloc)
                jindex, = np.where(diff < sliceeps)
                jindex = jindex[np.logical_and(jindex >= p.li[1], jindex < p.ui[1])]
                act_slice = x_i[1][jindex.min()]
                jindex  = slice(jindex.min(),jindex.max()+1)
            elif slicedir == 3 and kone > 0:
                diff = np.abs(x_i[2] - sliceloc)
                kindex, = np.where(diff < sliceeps)
                kindex = kindex[np.logical_and(kindex >= p.li[2], kindex < p.ui[2])]
                if len(kindex) == 0: continue # skip patch
                act_slice = x_i[2][kindex.min()]
                kindex  = slice(kindex.min(),kindex.max()+1)

            # establish integer offset for uniform plotting of vectors (across multiple patches)
            fvoffset = np.ceil((p.llc - global_llc)/p.dx)
            ivoffset = [(ivfreql - (int(o)) % ivfreql) % ivfreql for o in fvoffset]
            print "id: ", p.id, "lvl: ", p.level," VEC FREQ: ", ivfreql, "; ivoffset: ", ivoffset
            viindex = slice(iindex.start+ivoffset[0],min(iindex.stop,p.ui[0]),ivfreql)
            vjindex = slice(iindex.start+ivoffset[1],min(jindex.stop,p.ui[1]),ivfreql)
            vkindex = slice(kindex.start+ivoffset[2],min(kindex.stop,p.ui[2]),ivfreql)
            if slicedir == 1 or ione == 0:
                viindex = iindex
            elif slicedir == 2 or jone == 0:
                vjindex = jindex
            elif slicedir == 3 or kone == 0:
                vkindex = kindex

            times.append(p.time)

            if pltvar[0] == 'v1':
                var = p.data[iindex,jindex,kindex,p.varidx['p1']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
            elif pltvar[0] == 'v2':
                var = p.data[iindex,jindex,kindex,p.varidx['p2']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                var = var / gf['h2c'][iindex,np.newaxis,np.newaxis]
                if p.mesh_type == 'spherical':
                    var = var / x_i[0][iindex,np.newaxis,np.newaxis]
            elif pltvar[0] == 'v3':
                var = p.data[iindex,jindex,kindex,p.varidx['p3']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d']]
                var = ( var / gf['h31c'][iindex,np.newaxis,np.newaxis] / gf['h32c'][np.newaxis,jindex,np.newaxis] )
            elif pltvar[0] == 'divb':
                var, sumd, normmax = divergence(p, (p.varidx['b1'],p.varidx['b2'],p.varidx['b3']), gf)
                var = var[iindex,jindex,kindex]
            else:
                var = p.data[iindex,jindex,kindex,p.varidx[pltvar[0]]]
            dmass = p.data[p.li[0]-1:p.ui[0],p.li[1]-1:p.ui[1],p.li[2]-1:p.ui[2],p.varidx['d']] * p.ds.prod()
            mass = np.sum(dmass)

            # record global max. and min. of data
            qmin = min(qmin, np.nanmin(var))
            qmax = max(qmax, np.nanmax(var))

            # mask data for cutting
            if cut_mask.ndim == 0:
                var = np.ma.array(var,mask=cut_mask)
            else:
                var = np.ma.array(var,mask=cut_mask[iindex,jindex,kindex])
            cpltvar = set_variable_name(pltvar[0],passthrough=False)
            if vunits is not None:
                cpltvar = cpltvar + " [" + vunits[0] + "]"
            print cpltvar, ": ", var.min(), var.max(), var.shape, "; M = ", mass
            print iindex,jindex,kindex
            print viindex,vjindex,vkindex
            # Note that these values, in general, will differ between patches.
            if cmin[0] == 0.0 and cmax[0] == 0.0:
                cmin[0] = qmin
                cmax[0] = qmax

            if vscaling[0] is not None:
                var = var * vscaling[0]
                if m1 == 0 and m2 == 0:
                    cmin[0] = cmin[0] * vscaling[0]
                    cmax[0] = cmax[0] * vscaling[0]
            if clog[0] == 1:
                var = np.log10(var)
                cpltvar = r'log ' + cpltvar
            if m1 == 0 and m2 == 0 and clog[0] == 1:
                if cmin[0] != 0.0: cmin[0] = np.log10(cmin[0])
                if cmax[0] != 0.0: cmax[0] = np.log10(cmax[0])

            levels = np.linspace(cmin[0], cmax[0], num=50)
            conf = cax.contourf(xp[iindex,jindex,kindex].squeeze(),yp[iindex,jindex,kindex].squeeze(),
                                var.squeeze(),levels,cmap=ColMap,extend='both',zorder=-min(p.dx))

            if vecvar is not None and vecvar != '':
                if (viindex.step < viindex.stop) and (vjindex.step < vjindex.stop) and (vkindex.step < vkindex.stop):
                    if vecvar[0][0] == 'v':
                        vvar1 = p.data[viindex,vjindex,vkindex,p.varidx['p1']] / p.data[viindex,vjindex,vkindex,p.varidx['d']]
                        vvar2 = p.data[viindex,vjindex,vkindex,p.varidx['p2']] / p.data[viindex,vjindex,vkindex,p.varidx['d']]
                    elif vecvar[0][0] == 'p':
                        vvar1 = p.data[viindex,vjindex,vkindex,p.varidx[vecvar[0]]]
                        vvar2 = p.data[viindex,vjindex,vkindex,p.varidx[vecvar[1]]]
                    elif vecvar[0][0] == 'b':
                        vvar1 = p.data[viindex,vjindex,vkindex,p.varidx[vecvar[0]]]
                        vvar2 = p.data[viindex,vjindex,vkindex,p.varidx[vecvar[1]]]
                    vmagmax = (np.sqrt(vvar1**2 + vvar2**2)).max()
                    vmagmax_global = max(vmagmax_global, vmagmax)
                    print "Max. vector magnitude: ",vmagmax
                    quiv = cax.quiver(xp[viindex,vjindex,vkindex].squeeze(),
                                    yp[viindex,vjindex,vkindex].squeeze(),
                                    vvar1.squeeze(),vvar2.squeeze(),scale_units='width',scale=vpscale,
                                    color='k',zorder=-min(p.dx)*0.25,units='dots',width=1.0)

            if (showgrids or pltvar[0] == 'gs'):
                if xscaling is not None:
                    edges = patch_grid_edges(p,slicedir,scaling=xscaling,truecoords=truecoords)
                else:
                    edges = patch_grid_edges(p,slicedir,truecoords=truecoords)
                cax.add_patch(edges)
                #for c in conf.collections:
                #    c.set_clip_path(edges)

            # annotate patches with their ID (if desired)
            if showgridno:
                offset = [2*x for x in p.nghost]
                offset[slicedir-1] = 0
                gn_xpos = xp[iindex,jindex,kindex][tuple(offset)]
                gn_ypos = yp[iindex,jindex,kindex][tuple(offset)]
                cax.annotate(str(p.id),(gn_xpos,gn_ypos),color='k',weight='extra bold',
                            fontsize=13,zorder=min(p.dx))

            m2 += 1

        cbar = plt.colorbar(conf, ax=cax, fraction=0.03, pad=0.03)
        cbar.ax.tick_params(axis='both',which='major',labelsize=16)
        cbar.ax.annotate(cpltvar,xy=(-0.10,-0.06),fontsize=14,xycoords='axes fraction')
        if vecvar is not None and vecvar != '':
            quivlength = 1.0
            exponent = int(np.log10(quivlength*vecscaling))
            prefactor = quivlength*vecscaling / 10**exponent
            fmt_string = "{0:4.3g}".format(prefactor) + r"$\times 10^" + "{0:0d}$ ".format(exponent) + vecunits
            quivkey = cax.quiverkey(quiv,1.11,0.02,quivlength,fmt_string,
                                    coordinates='axes',labelpos='S',fontproperties={'size':12})

        avg_time = np.mean(times)
        if not suppress_title:
            cfig.suptitle(r"$\langle t \rangle$ = {0:.3e}; slice = {1:.3e}".format(avg_time,act_slice),fontsize=16)
        cfig.tight_layout(pad=0.7)

        if output_png:
            currfile = "p{0:04d}.png".format(m1)
            print currfile
            cfig.savefig(currfile,dpi=100)
            plt.close(cfig)
        else:
            plt.draw()

        m1 += 1
        print "qmin, qmax: ", qmin, qmax
        print "vmagmax: ", vmagmax_global
        print '================'
