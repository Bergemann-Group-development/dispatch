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

def patch_edges(p, slicedir, pcolr, truecoords=True):
    """
    Using the corners of the patch, create a matplotlib patch that draws
    the grid edges.

    """

    from matplotlib.patches import Rectangle, Wedge

    assert isinstance(p, dd.Patch), "Input must be an instance of a DISPATCH `Patch`!"

    if truecoords:
        if slicedir == 1:
            llx, urx = p.active_corners[0][1], p.active_corners[1][1]
            lly, ury = p.active_corners[0][2], p.active_corners[1][2]
        elif slicedir == 2:
            llx, urx = p.active_corners[0][0], p.active_corners[1][0]
            lly, ury = p.active_corners[0][2], p.active_corners[1][2]
        elif slicedir == 3:
            llx, urx = p.active_corners[0][0], p.active_corners[1][0]
            lly, ury = p.active_corners[0][1], p.active_corners[1][1]
    else:
        if slicedir == 1:
            llx, urx = p.llc[1], p.llc[1] + p.size[1]
            lly, ury = p.llc[2], p.llc[2] + p.size[2]
        elif slicedir == 2:
            llx, urx = p.llc[0], p.llc[0] + p.size[0]
            lly, ury = p.llc[2], p.llc[2] + p.size[2]
        elif slicedir == 3:
            llx, urx = p.llc[0], p.llc[0] + p.size[0]
            lly, ury = p.llc[1], p.llc[1] + p.size[1]

    if pp.id == 1:
        llx *= 1.08
        lly *= 1.08
        urx *= 1.08
        ury *= 1.08
    if p.mesh_type == dg.mesh_types[1] or truecoords: # Cartesian
        edges = Rectangle((llx,lly),width=(urx-llx),height=(ury-lly),fc=pcolr,lw=1.5,
                          ec='0.2',linestyle='solid',zorder=-min(p.dx))

    return edges

if __name__ == "__main__":

    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)
    plt.interactive(True)

    data = './data/input_nested'
    iters = [0]
    pxsize = [-0.57, 0.57]
    pysize = [-0.57, 0.57]
    showgridno = False
    ColMap = cm.RdBu_r
    # for best results, skip patches which are entirely overlapped by finer patches
    skipped = [6,15]

    pcolrs=['forestgreen','royalblue','salmon','orchid']

    mi = 0
    for i in iters:
        if not os.path.exists(data):
            raise Exception("data dir not found!")
        patches = du.patches(data=data,iout=i)
        if len(patches) == 0: raise Exception("no patches found!")
        for pp in patches:
            pp.cache()
            pp.stripit()

        cfig = plt.figure("{0}".format(i),facecolor='w',figsize=(4.9,5.0))
        cfig.clf()
        gs = gridspec.GridSpec(1,1)

        axs = []
        for j in xrange(1):
            if j == 0:
                cax = cfig.add_subplot(gs[0,j],aspect=1.0)
            else:
                cax = cfig.add_subplot(gs[0,j],sharex=axs[0],sharey=axs[0])
            axs.append(cax)
            cax.axis(xmin=pxsize[0],xmax=pxsize[1],ymin=pysize[0],ymax=pysize[1])
            cax.xaxis.set_ticks_position('none')
            cax.yaxis.set_ticks_position('none')
            cax.tick_params(labelbottom=False)
            cax.tick_params(labelleft=False)
            for loc,spine in cax.spines.items():
                if loc in ['right','top','left','bottom']:
                    spine.set_color('none')
            # ensure the axes remain visible by putting them at the "top"
            cax.set_zorder(99.0)

            maxlevel = 0
            for pp in patches:
                try:
                    skipped.index(pp.id)
                    print "Skipping ",pp.id," (user)"
                    continue
                except:
                    pass

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

                print "Data:", pp.filename, "; step #: ",pp.istep

                edges = patch_edges(pp,3,pcolrs[pp.level-1],truecoords=True)
                cax.add_patch(edges)
                # annotate patches with their ID (if desired)
                if showgridno:
                    offset = tuple([2*x for x in pp.nghost])
                    gn_xpos = xp[iindex,jindex,kindex][offset]
                    gn_ypos = yp[iindex,jindex,kindex][offset]
                    cax.annotate(str(pp.id),(gn_xpos,gn_ypos),color='k',weight='extra bold',
                                 fontsize=14,zorder=10.0)

        cax.annotate(1,(-0.53,-0.53),color='k',weight='extra bold',fontsize=15,zorder=10)
        cax.annotate(2,(-0.45,-0.45),color='k',weight='extra bold',fontsize=15,zorder=10)
        cax.annotate(3,(-0.14,-0.14),color='k',weight='extra bold',fontsize=15,zorder=10)
        cax.annotate(4,(-0.05,-0.05),color='k',weight='extra bold',fontsize=15,zorder=10)
        cfig.subplots_adjust(bottom=0.095,left=0.0491,right=0.951,top=0.97)
        plt.draw()
        print "="*20

        # increment counter (`iters` loop)
        mi += 1

