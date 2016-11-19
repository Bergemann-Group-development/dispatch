import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from dispatch_data import *
import itertools
#import ipdb

def dump_min_max(p):
    """
    Print the min and max of each variable in `p.data` (convenience function).

    """

    print "Time: {0}".format(p.time)
    for q,v in zip(p.varidx.iterkeys(),p.varidx.itervalues()):
        print "{0}: min = {1}; max = {2}".format(q,p.data[:,:,:,v].min(),p.data[:,:,:,v].max())

def plotting_grid(root, p, x_i, slicedir, sliceloc, sliceeps=1.0e-2, truecoords=False,
                  rootcoords=False, phi_rotate=0.0, theta_rotate=0.0):
    """
    Given slice information, return the appropriate coord. vectors for plotting,
    including the conversion to Cartesian coords. (required for plotting with
    matplotlib).
    
    If the `truecoords` keyword is set True, then the conversion to Cartesian
    coords. is skipped.
    
    Note that `slicedir` and `sliceloc` depend on `truecoords`: if `truecoords`
    is False, they should be in Cartesian coords., otherwise the should be in
    native coord. system of the patch.
    
    The following curvilinear coord. systems are known:
      (x,y,z), (Z,R,Phi), (r,theta,phi)
    The `phi_rotate` and `theta_rotate` arguments allow for patch coords. which
    are rotated relative to the convention of theta = phi = 0. They are silently
    passed to the coord. transformation functions.

    """

    assert (slicedir > 1) or (slicedir < 3), "slice direction must be 1, 2, or 3."
    assert not (truecoords and rootcoords), "truecoords and rootcoords parameters are mutually exclusive."

    if p.mesh_type == 'Cartesian' or truecoords:
        x, y, z = np.meshgrid(x_i[0], x_i[1], x_i[2], indexing='ij')
    elif p.mesh_type == 'cylindrical':
        x, y, z = zrp_to_xyz(x_i[0], x_i[1], x_i[2], phi_rotate=phi_rotate)
    elif p.mesh_type == 'spherical':
        x, y, z = rtp_to_xyz(x_i[0], x_i[1], x_i[2], phi_rotate=phi_rotate, theta_rotate=theta_rotate)

    # adjust coords. using the patch position (which is in Cartesian coords.).
    #x += p.pos[0]
    #y += p.pos[1]
    #z += p.pos[2]

    # assign coords. actually needed for plotting; cyclical.
    if slicedir == 1:
        xp, yp = y,z
    elif slicedir == 2:
        xp, yp = z,x
    elif slicedir == 3:
        xp, yp = x,y

    # account for symmetry (i.e., less than 3-D).
    sliceloc_act = sliceloc
    if p.mesh_type == 'Cartesian':
        if slicedir == 1 and len(x_i[0]) == 1: sliceloc_act = x_i[0]
        if slicedir == 2 and len(x_i[1]) == 1: sliceloc_act = x_i[1]
        if slicedir == 3 and len(x_i[2]) == 1: sliceloc_act = x_i[2]
    if p.mesh_type == 'cylindrical':
        if slicedir == 3 and len(x_i[0]) == 1: sliceloc_act = x_i[0]
    if p.mesh_type == 'spherical':
        if slicedir == 3 and len(x_i[1]) == 1: sliceloc_act = np.cos(x_i[1])

    # finally, determine the (nearest) slice location for the current patch.
    if slicedir == 1:
        diff = x - sliceloc_act
        sliceindices = np.where(abs(diff) < sliceeps)[0]
    elif slicedir == 2:
        diff = y - sliceloc_act
        sliceindices = np.where(abs(diff) < sliceeps)[1]
    elif slicedir == 3:
        diff = z - sliceloc_act
        sliceindices = np.where(abs(diff) < sliceeps)[2]
    sliceindices = set(sliceindices)
    # create mask
    m = np.ma.masked_outside(diff,-sliceeps,sliceeps)
    xpcut = np.ma.array(xp,mask=m.mask)
    ypcut = np.ma.array(yp,mask=m.mask)
    
    return xpcut, ypcut, m.mask, sliceindices

def patch_grid_edges(p):
    """
    Using the corners of the patch, create a matplotlib patch that draws
    the grid edges.

    """

    from matplotlib.patches import Rectangle, Wedge

    assert isinstance(p, Patch), "Input must be an instance of a DISPATCH `Patch`!"

    width = p.active_corners[1][0] - p.active_corners[0][0]
    height = p.active_corners[1][1] - p.active_corners[0][1]

    if p.mesh_type == mesh_types[1]: # Cartesian
        edges = Rectangle(p.active_corners[0][0:2],width,height,fc='none',lw=1.0,
                          ec='k',linestyle='dashed',zorder=0.01/min(p.dx)+0.1)
    elif p.mesh_type == mesh_types[2]: # spherical
        thetamin = round(p.active_corners[0][2]*180.0/np.pi,3)
        thetamax = round(p.active_corners[1][2]*180.0/np.pi,3)
        llx, urx = p.active_corners[0][0], p.active_corners[1][0]
        edges = Wedge((p.pos[0],p.pos[1]),urx,thetamin,thetamax,width=(urx-llx),fc='none',lw=1.0,
                      ec='k',linestyle='dotted',zorder=0.01/min(p.dx)+0.1)
    elif p.mesh_type == mesh_types[3]: # cylindrical
        phimin = round(p.active_corners[0][2]*180.0/np.pi,3)
        phimax = round(p.active_corners[1][2]*180.0/np.pi,3)
        llx, urx = p.active_corners[0][1], p.active_corners[1][1]
        edges = Wedge((p.pos[0],p.pos[1]),urx,phimin,phimax,width=(urx-llx),fc='none',lw=1.0,
                      ec='k',linestyle='dotted',zorder=0.01/min(p.dx)+0.1)

    return edges

def dispatch_grid(p):
    """
    Given grid information for a patch, calculate the three *cell-centred*
    coordinate arrays in the coord. system of the patch.

    """

    if p.n != p.nwrite:
        x1 = p.corners[0][0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
        x2 = p.corners[0][1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
        x3 = p.corners[0][2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]
    else:
        x1 = p.active_corners[0][0] + np.arange(p.nwrite[0]) * p.dx[0] + 0.5 * p.dx[0]
        x2 = p.active_corners[0][1] + np.arange(p.nwrite[1]) * p.dx[1] + 0.5 * p.dx[1]
        x3 = p.active_corners[0][2] + np.arange(p.nwrite[2]) * p.dx[2] + 0.5 * p.dx[2]

    return (x1, x2, x3)

if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    gamma = 5.0/3.0
    datadir = './data/'
    prefix = ''
    ptimes = [2.0e-3]
    t_eps = 0.03
    pxsize = [0.0,1.0]
    pysize = [0.0,1.0]
    #pltvar = ['d']
    pltvar = ['J']
    cmin   = [0.0]
    cmax   = [5.0]
    clog   = [  0]
    slicedir = 3  # given with respect to the root grid (typically patch id = 1)
    sliceloc = 0.501 # coord. of slice to be taken
    showgrids = True
    showgridno = False
    output_png = False
    truecoords = False
    ThisPatch = StaggerMHDPatch
    ColMap = cm.RdBu_r
    verbose = False
    skipped = []

    m1 = 0 # counter 1
    root = None
    for pt in ptimes:

        cfig = plt.figure("t = {0:.3e}".format(pt),facecolor='w',figsize=(7.0,6.0))
        cfig.clf()
        cax = cfig.add_subplot(1,1,1)
        cax.axis(xmin=pxsize[0],xmax=pxsize[1],ymin=pysize[0],ymax=pysize[1])
        #cax.axes.set_aspect('equal')
        cax.set_xlabel(r'$x$')
        cax.set_ylabel(r'$y$')
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
            if (abs(ptemp.time-pt)/(pt+1.0e-99) <= t_eps): patches.append(ptemp)
            if (verbose): print '----'
        if len(patches) == 0:
            print "No patches found?! t = ",pt
            continue

        qmin, qmax = 1.0e99, -1.0e99
        m2 = 0 # counter 2
        for p in patches:
            if p.task_type == 'immersed_boundary_t':
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

            print p.filename, p.time
            times.append(p.time)

            # hack in the mean intensity as a plottable variable
            p.varidx['J'] = p.nvar-1

            x_i = dispatch_grid(p)
            xp, yp, cut_mask, slice_index = plotting_grid(root, p, x_i, slicedir, sliceloc,
                                               truecoords=truecoords,sliceeps=0.006)
            if (np.ma.count_masked(xp) == np.product(p.nwrite) or
                np.ma.count_masked(yp) == np.product(p.nwrite)):
                   print "Skipping ID {0}..".format(p.id)
                   continue
            if len(slice_index) > 1:
                print "WARNING!: More than one index matches the desired slice!"
                print slice_index

            # ghost zones *can* now be included (or not)
            ione = (0,1)[p.n[0] > 1]
            jone = (0,1)[p.n[1] > 1]
            kone = (0,1)[p.n[2] > 1]
            if p.n != p.nwrite:
                iindex = slice(p.nghost[0] - ione,p.n[0] + p.nghost[0] + ione)
                jindex = slice(p.nghost[1] - jone,p.n[1] + p.nghost[1] + jone)
                kindex = slice(p.nghost[2] - kone,p.n[2] + p.nghost[2] + kone)
            else:
                iindex = slice(0,p.n[0])
                jindex = slice(0,p.n[1])
                kindex = slice(0,p.n[2])

            if pltvar[0] == 'v1':
                var = p.data[iindex,jindex,kindex,p.varidx['s1']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d ']]
            elif pltvar[0] == 'v2':
                var = p.data[iindex,jindex,kindex,p.varidx['s2']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d ']]
                if p.mesh_type == 'spherical':
                    var = var / x_i[0][iindex,np.newaxis,np.newaxis]
            elif pltvar[0] == 'v3':
                var = p.data[iindex,jindex,kindex,p.varidx['s3']]
                var = var / p.data[iindex,jindex,kindex,p.varidx['d ']]
                if p.mesh_type == 'cylindrical':
                    var = var / x_i[1][np.newaxis,jindex,np.newaxis]
                if p.mesh_type == 'spherical':
                    var = ( var / x_i[0][iindex,np.newaxis,np.newaxis]
                          / np.sin(x_i[1][np.newaxis,jindex,np.newaxis]) )
            else:
                var = p.data[iindex,jindex,kindex,p.varidx[pltvar[0]]]

            # record global max. and min. of data
            qmin = min(qmin, np.nanmin(var))
            qmax = max(qmax, np.nanmax(var))

            # mask data for cutting
            var = np.ma.array(var,mask=cut_mask[iindex,jindex,kindex])
            cpltvar = pltvar[0]
            print cpltvar, ": ", var.min(), var.max(), var.shape
            if clog[0] == 1:
                var = np.log10(var)
                cpltvar = r'log ' + cpltvar
            if m1 == 0 and m2 == 0 and clog[0] == 1:
                if cmin[0] != 0.0: cmin[0] = np.log10(cmin[0])
                if cmax[0] != 0.0: cmax[0] = np.log10(cmax[0])

            levels = np.linspace(cmin[0], cmax[0], num=50)
            

            if slicedir == 1:
                iindex = slice_index.pop()
                vvv = var[iindex-1,:,:]
            elif slicedir == 2:
                jindex = slice_index.pop()
                vvv = var[:,jindex-1,:]
            elif slicedir == 3:
                kindex = slice_index.pop()
                vvv = var[:,:,kindex-1]
            xpp = xp[iindex,jindex,kindex]
            ypp = yp[iindex,jindex,kindex]
            #ipdb.set_trace()
            conf = cax.contourf(xpp,ypp,vvv
                                ,levels,cmap=ColMap,extend='both',
                                zorder=-min(p.dx))

            if (showgrids or pltvar[0] == 'gs') and p.id != 1:
                edges = patch_grid_edges(p)
                cax.add_patch(edges)
                #for c in conf.collections:
                #    c.set_clip_path(edges)

            # annotate patches with their ID (if desired)
            if p.id != 1 and showgridno:
                offset = [3*x for x in p.nghost]
                offset.pop(slicedir-1)
                gn_xpos = xpp[offset]
                gn_ypos = ypp[offset]
                cax.annotate(str(p.id),(gn_xpos,gn_ypos),color='k',weight='extra bold',
                            fontsize=13*root.level/p.level,zorder=min(p.dx))

            m2 += 1

        cbar = plt.colorbar(conf, ax=cax, fraction=0.03, pad=0.03)
        cbar.ax.tick_params(axis='both',which='major',labelsize=16)
        cbar.ax.annotate(cpltvar,xy=(-0.10,-0.06),fontsize=14,xycoords='axes fraction')

        avg_time = np.mean(times)
        cfig.suptitle(r"$\langle t \rangle$ = {0:.3e}".format(avg_time),fontsize=16)
        cfig.tight_layout()

        if output_png:
            currfile = "p{0:04d}.png".format(iout)
            print currfile
            cfig.savefig(currfile,dpi=100)
            plt.close(cfig)
        else:
            plt.draw()

        m1 += 1
        print "qmin, qmax: ", qmin, qmax
        print '================'
