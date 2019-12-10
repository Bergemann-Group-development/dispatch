import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
import itertools

if __name__ == "__main__":

    plt.interactive(True)
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=18)

    # single resolution, single patch
    single_ncells   = [     16.,     32.,     64.,    128.,    256.,    512.,   1024.,   2048.]
    single_Linf_max = [1.649e-2,5.917e-3,1.541e-3,3.873e-4,9.651e-5,2.415e-5,6.032e-6,1.599e-6]

    # multi-patch, single resolution; each patch is 16**2 or 32**2 cells
    # 3^2 patches
    multi_n_per_patch3  = [       8,      16,      32,      64,     128,     256]
    multi_equiv_ncells3 = [      24,      48,      96,     192,     384,     768]
    multi_Linf_max3     = [1.635e-2,4.116e-3,9.898e-4,2.370e-4,5.905e-5,1.511e-5]
    
    # same, but 9^2 patches
    mulit_n_per_patch9  = [       4,       8,      16,      32,       64,     128]
    multi_equiv_ncells9 = [      36,      72,     144,     288,      576,    1152]
    multi_Linf_max9     = [1.542e-2,3.612e-3,8.595e-4,2.157e-4, 5.574e-5,1.482e-5]

    # nested multi-patch
    # a refinement ratio of 3 is used
    # 2 levels
    #nested_2l_equiv_ncells = [      48,      96,     144,     288]
    #nested_2l_Linf_max     = [6.867e-3,2.232e-3,7.219e-4,1.630e-4]
    # 3 levels
    nested_3l_n_per_patch  = [       5,       9,      17,      33,      65,     129]
    nested_3l_equiv_ncells = [      45,      81,     153,     297,     585,    1161]
    nested_3l_Linf_max     = [9.234e-3,1.925e-3,5.262e-4,1.300e-4,3.275e-5,8.444e-6]
    # 4 levels
    nested_4l_n_per_patch  = [       9,      17,      33,      65]
    nested_4l_equiv_ncells = [     243,     459,     891,    1755]
    nested_4l_Linf_max     = [2.448e-4,5.470e-5,1.376e-5,3.429e-6]

    # mock data for convergence-order lines
    xcells = np.logspace(4.0,10.2,base=2.0)
    offset = 2.7

    fig1 = plt.figure("L-infinity static potential",facecolor='w',figsize=(7.5,5.0))
    fig1.clf()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.set_ylabel(r"$L_\infty = \Vert\Phi - \Phi_\mathrm{e}\Vert_\infty / \Vert\Phi_\mathrm{e}\Vert_\infty$")
    ax1.set_xlabel(r"$N$, effective \# of cells per direction")
    ax1.axis(xmin=1e1,xmax=3e3,ymin=1e-6,ymax=3e-2)

    ms = 9
    ax1.loglog(single_ncells,single_Linf_max,'o',ms=ms,mew=0.0,label='single patch',color='forestgreen')
    ax1.loglog(multi_equiv_ncells3,multi_Linf_max3,'s',ms=ms,mew=0.0,label='$3^2$ patches',color='royalblue')
    ax1.loglog(multi_equiv_ncells9,multi_Linf_max9,'s',ms=ms,mew=0.0,label='$9^2$ patches',color='c')
    ax1.loglog(nested_3l_equiv_ncells,nested_3l_Linf_max,'^',ms=ms+1.0,mew=0.0,label='nested, $l_\mathrm{max} = 3$',color='salmon')
    ax1.loglog(nested_4l_equiv_ncells,nested_4l_Linf_max,'^',ms=ms+1.0,mew=0.0,label='nested, $l_\mathrm{max} = 4$',color='orchid')
    ax1.loglog(xcells,offset*xcells**(-2),'k--',lw=1.0)
    ax1.annotate(xy=(1e2,9e-5),s=r"$\propto 1/N^2$",fontsize=16,rotation=-43.5)

    leg = ax1.legend(loc='upper right',frameon=True,numpoints=1,handletextpad=0.0,fontsize=17,borderaxespad=0.0)
    lframe = leg.get_frame()
    lframe.set_facecolor('w')
    lframe.set_edgecolor('none')

    # aesthetics
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    for loc,spine in ax1.spines.items():
        if loc in ['right','top']:
            spine.set_color('none')

    ax1.grid(True)

    fig1.tight_layout()
    fig1.subplots_adjust(left=0.13,right=0.96,bottom=0.12,top=0.95)
    plt.draw()
