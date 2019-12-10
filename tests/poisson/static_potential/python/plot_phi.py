import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, patches as mplpatches, gridspec
import dispatch
import itertools
import fnmatch
import os

def patch_grid_edges(p, slicedir, scale_factor=1.0, truecoords=True):
    """
    Using the corners of the patch, create a matplotlib patch that draws
    the grid edges.

    """

    from matplotlib.patches import Rectangle, Wedge

    lo_corners = []
    hi_corners = []
    origin = p.llc_cart
    for o,s in zip(origin,p.size):
        lo_corners.append(o)
        hi_corners.append(o + s)
    active_corners = lo_corners, hi_corners

    if truecoords:
        if slicedir == 1:
            llx, urx = active_corners[0][1], active_corners[1][1]
            lly, ury = active_corners[0][2], active_corners[1][2]
        elif slicedir == 2:
            llx, urx = active_corners[0][0], active_corners[1][0]
            lly, ury = active_corners[0][2], active_corners[1][2]
        elif slicedir == 3:
            llx, urx = active_corners[0][0], active_corners[1][0]
            lly, ury = active_corners[0][1], active_corners[1][1]
    else:
        if slicedir == 1:
            llx, urx = p.llc_cart[1], p.llc_cart[1] + p.size[1]
            lly, ury = p.llc_cart[2], p.llc_cart[2] + p.size[2]
        elif slicedir == 2:
            llx, urx = p.llc_cart[0], p.llc_cart[0] + p.size[0]
            lly, ury = p.llc_cart[2], p.llc_cart[2] + p.size[2]
        elif slicedir == 3:
            llx, urx = p.llc_cart[0], p.llc_cart[0] + p.size[0]
            lly, ury = p.llc_cart[1], p.llc_cart[1] + p.size[1]

    if p.mesh_type == 'Cartesian' or truecoords: # Cartesian or Cartesian-like
        llx *= scale_factor
        urx *= scale_factor
        if p.mesh_type == 'Cartesian':
            lly *= scale_factor
            ury *= scale_factor
        edges = Rectangle((llx,lly),width=(urx-llx),height=(ury-lly),fc='none',lw=1.0,
                          ec='0.7',linestyle='dashed',zorder=0.01/min(p.ds)+0.1)
    elif p.mesh_type == 'spherical':
        # note that matplotlib sets theta=0 to what we typically consider to be theta=pi/2.
        thetamin = round(active_corners[0][1]*180.0/np.pi-90.0,3)
        thetamax = round(active_corners[1][1]*180.0/np.pi-90.0,3)
        llx, urx = active_corners[0][0], active_corners[1][0]
        if slicedir == 3 or slicedir == 2:
            llx *= scale_factor
            urx *= scale_factor
        edges = Wedge((0.0,0.0),urx,thetamin,thetamax,width=(urx-llx),fc='none',
                      lw=1.0,ec='k',linestyle='dotted',zorder=0.01/min(p.ds)+0.1)
    elif p.mesh_type == 'cylindrical':
        llx *= scale_factor
        urx *= scale_factor
        if slicedir == 3:
            lly *= scale_factor
            ury *= scale_factor
        # note that matplotlib sets theta=0 to what we typically consider to be theta=pi/2.
        phimin = round(active_corners[0][2]*180.0/np.pi-90.0,3)
        phimax = round(active_corners[1][2]*180.0/np.pi-90.0,3)
        llx, urx = active_corners[0][1], active_corners[1][1]
        edges = Wedge((p.position[0]*scale_factor,p.position[1]*scale_factor),urx,phimin,phimax,width=(urx-llx),fc='none',
                      lw=1.0,ec='k',linestyle='dotted',zorder=0.01/min(p.ds)+0.1)

    return edges

def phi_exact(x,y,z,r0=1.0):
    """
    \Phi_e(r) = ln( (r/r_0)^2 + 1)

    """

    if x.ndim == 1 or y.ndim == 1 or z.ndim == 1:
        xx, yy, zz = np.meshgrid(x,y,z,indexing='ij')
    else:
        xx, yy, zz = x, y, z
    r = np.sqrt(xx**2 + yy**2)
    phi_e = np.log((r / r0)**2 + 1.0) + 1.0

    return phi_e

def rho_exact(x,y,z,r0=1.0):
    """
    \rho_e(r) = 4 * r_0^2 / (r^2 + r_0^2)^2

    """

    if x.ndim == 1 or y.ndim == 1 or z.ndim == 1:
        xx, yy, zz = np.meshgrid(x,y,z,indexing='ij')
    else:
        xx, yy = x, y, z
    r = np.sqrt(xx**2 + yy**2)
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

    data = '../data'
    iters = [0,10,50]
    #iters = [0]
    pxsize = [-0.5, 0.5]
    pysize = [-0.5, 0.5]
    pltvar = ['rdphi','phie', 'phi', 'rdd']
    cmin   = [   1e-4,   1.0,   0.0,  1e-6]
    cmax   = [    1.0,    5.,    5.,   0.5]
    clog   = [      1,     1,     1,     1]
    showgrids = True
    showgridno = True
    ColMap = cm.RdBu_r
    # for best results, skip patches which are entirely overlapped by finer patches
    skipped = []
    r0 = 0.1

    L2_min, L2_max, L2_mean = [], [], []
    Linf_min, Linf_max, Linf_mean = [], [], []
    istep_min, istep_max = [], []
    mi = 0
    for i in iters:
        if not os.path.exists(data):
            raise Exception("data dir not found!")
        snap = dispatch.snapshot(data=data,iout=i)
        patches = snap.patches
        if len(patches) == 0: raise Exception("no patches found!")

        cfig = plt.figure("{0}".format(i),facecolor='w',figsize=(5.0*len(pltvar),4.0))
        cfig.clf()
        gs = gridspec.GridSpec(1,len(pltvar))
        
        axs = []
        istep_p_min, istep_p_max = 1000000, -1000000
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

            spltvar = pltvar[j]
            if clog[j] == 1:
                spltvar = r'log ' + spltvar
                if cmin[j] != 0.0 and mi == 0: cmin[j] = np.log10(cmin[j])
                if cmax[j] != 0.0 and mi == 0: cmax[j] = np.log10(cmax[j])
            levels = np.linspace(cmin[j],cmax[j],num=20)
            print cmin[j],cmax[j]

            for pp in patches:
                try:
                    skipped.index(pp.id)
                    print "Skipping ",pp.id," (user)"
                    continue
                except:
                    pass
                # ghost zones *can* now be included (or not)
                ione = (0,1)[pp.n[0] > 1]
                jone = (0,1)[pp.n[1] > 1]
                kone = (0,1)[pp.n[2] > 1]
                if np.any(pp.n != pp.ncell):
                    iindex = slice(pp.ng[0] - ione,pp.n[0] + pp.ng[0] + ione)
                    jindex = slice(pp.ng[1] - jone,pp.n[1] + pp.ng[1] + jone)
                    kindex = slice(pp.ng[2] - kone,pp.n[2] + pp.ng[2] + kone)
                else:
                    iindex = slice(0,pp.n[0])
                    jindex = slice(0,pp.n[1])
                    kindex = slice(0,pp.n[2])

                x_i = (pp.x, pp.y, pp.z)
                xp = pp.x[iindex]
                yp = pp.y[jindex]
                zp = pp.z[kindex]

                print "Data:", pp.filename, "; step #: ",pp.istep
                if pltvar[j] == 'rdd':
                    var = pp.var('d',all=True)[iindex,jindex,kindex]
                    q0 = rho_exact(pp.x,pp.y,pp.z,r0=r0)
                    var = relative_difference(q0[iindex,jindex,kindex].squeeze(),var.squeeze())
                elif pltvar[j] == 'rdphi':
                    var = pp.var('phi',all=True)[iindex,jindex,kindex]
                    q0 = phi_exact(pp.x,pp.y,pp.z,r0=r0)
                    var = relative_difference(q0[iindex,jindex,kindex].squeeze(),var.squeeze())
                elif pltvar[j] == 'dphi':
                    var = pp.var('phi',all=True)[iindex,jindex,kindex]
                    q0 = phi_exact(pp.x,pp.y,pp.z,r0=r0)
                    var = var.squeeze() - q0[iindex,jindex,kindex].squeeze()
                elif pltvar[j] == 'de':
                    var = rho_exact(pp.x,pp.y,pp.z,r0=r0)
                    var = var[iindex,jindex,kindex]
                elif pltvar[j] == 'phie':
                    var = phi_exact(pp.x,pp.y,pp.z,r0=r0)
                    var = var[iindex,jindex,kindex]
                elif pltvar[j] == 'dd':
                    var = pp.var('d',all=True)[iindex,jindex,kindex]
                    q0 = rho_exact(pp.x,pp.y,pp.z,r0=r0)
                    var = var.squeeze() - q0[iindex,jindex,kindex].squeeze()
                else:
                    var = pp.var(pltvar[j],all=True)[iindex,jindex,kindex]

                print spltvar, ": ", var.min(), var.max(), var.shape
                if j == 0:
                    istep_p_min = min(istep_p_min, pp.istep)
                    istep_p_max = max(istep_p_max, pp.istep)

                if clog[j] == 1:
                    var = np.log10(var)

                conf = cax.contourf(xp,yp,var.T.squeeze(),levels,cmap=ColMap,
                             extend='both',zorder=-min(pp.ds))

                if (showgrids or pltvar[0] == 'gs'):
                    edges = patch_grid_edges(pp,3,truecoords=True)
                    cax.add_patch(edges)
                # annotate patches with their ID (if desired)
                if showgridno:
                    offset = tuple([2*x for x in pp.ng])
                    gn_xpos = xp[offset[0]]
                    gn_ypos = yp[offset[1]]
                    cax.annotate(str(pp.id),(gn_xpos,gn_ypos),color='k',weight='extra bold',
                                 fontsize=13,zorder=min(pp.ds))

            cbar = plt.colorbar(conf, ax=cax, fraction=0.03, pad=0.03)
            cbar.ax.annotate(spltvar,xy=(-0.10,-0.06),fontsize=14,xycoords='axes fraction')

        istep_min.append(istep_p_min)
        istep_max.append(istep_p_max)

        cfig.tight_layout()
        plt.draw()
        print "="*20

        # print L2 and Linf norm of the result
        L2_iter, Linf_iter = [], []
        for pp in patches:
            try:
                skipped.index(pp.id)
                print "Skipping ",pp.id," (user)"
                continue
            except:
                pass

            # ghost zones *can* now be included (or not)
            ione = (0,1)[pp.n[0] > 1]
            jone = (0,1)[pp.n[1] > 1]
            kone = (0,1)[pp.n[2] > 1]
            if np.any(pp.n != pp.ncell):
                iindex = slice(pp.ng[0] - ione,pp.n[0] + pp.ng[0] + ione)
                jindex = slice(pp.ng[1] - jone,pp.n[1] + pp.ng[1] + jone)
                kindex = slice(pp.ng[2] - kone,pp.n[2] + pp.ng[2] + kone)
            else:
                iindex = slice(0,pp.n[0])
                jindex = slice(0,pp.n[1])
                kindex = slice(0,pp.n[2])

            x_i = (pp.x, pp.y, pp.z)
            xp = pp.x[iindex]
            yp = pp.y[jindex]
            zp = pp.z[kindex]

            phie = phi_exact(pp.x,pp.y,pp.z,r0=r0)
            dphi = pp.var('phi',all=True)[:,:,:] - phie
            L2d = np.linalg.norm(phie[iindex,jindex,kindex].flatten(),ord=2)
            L2n = np.linalg.norm(dphi[iindex,jindex,kindex].flatten(),ord=2)
            Linfd = np.linalg.norm(phie[iindex,jindex,kindex].flatten(),ord=np.inf)
            Linfn = np.linalg.norm(dphi[iindex,jindex,kindex].flatten(),ord=np.inf)
            L2_iter.append(L2n/L2d)
            Linf_iter.append(Linfn/Linfd)
            print "i = {0}, p = {1}, n = {2}, L2-norm = {3:.3e}, Linf-norm = {4:.3e}".format(i,pp.id,pp.n,L2n/L2d,Linfn/Linfd)

        # accumulate L2, Linf results
        L2_min.append(np.min(L2_iter))
        L2_max.append(np.max(L2_iter))
        L2_mean.append(np.mean(L2_iter))
        Linf_min.append(np.min(Linf_iter))
        Linf_max.append(np.max(Linf_iter))
        Linf_mean.append(np.mean(Linf_iter))

        # increment counter (`iters` loop)
        mi += 1

    for i in xrange(len(iters)):
        print "i = {0}, iout = {1}, min nsteps = {2}, max nsteps = {3}".format(i,iters[i],istep_min[i], istep_max[i])
        print "L2 min = {0:.3e}, L2 max = {1:.3e}, L2 mean = {2:.3e}".format(L2_min[i],L2_max[i],L2_mean[i])
        print "Linf min = {0:.3e}, Linf max = {1:.3e}, Linf mean = {2:.3e}".format(Linf_min[i],Linf_max[i],Linf_mean[i])
        print "-"*40
