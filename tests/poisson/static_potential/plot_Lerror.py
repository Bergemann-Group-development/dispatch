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
    multi_n_per_patch3  = [       8,      16,      32,      64,     128]
    multi_equiv_ncells3 = [      24,      48,      96,     192,     384]
    multi_Linf_max3     = [1.635e-2,4.116e-3,9.898e-4,3.695e-4,9.123e-5]
    
    # same, but 9^2 patches
    mulit_n_per_patch9  = [       8,      16,      32,       64,     128]
    multi_equiv_ncells9 = [      72,     144,     288,      576,    1152]
    multi_Linf_max9     = [4.156e-3,9.190e-4,6.132e-4, 1.663e-4,4.536e-5]

    # nested multi-patch
    # each patch is 32**2 cells, a refinement ratio of three is used
    # 2 levels
    nested_2l_equiv_ncells = [      36,      72,     144,     288,     216]
    nested_2l_Linf_max     = [4.739e-2,8.491e-3,7.219e-4,1.630e-4,1.641e-2]
    # 3 levels
    nested_3l_equiv_ncells = [     216,     432]
    nested_3l_Linf_max     = [1.641e-2,2.245e-3]

    # mock data for convergence-order lines
    xcells = np.logspace(3.6,11.5,base=2.0)
    offset = 2.5

    fig1 = plt.figure("L-infinity static potential",facecolor='w',figsize=(7.5,5.0))
    fig1.clf()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.set_ylabel(r"$L_\infty = \Vert\Phi - \Phi_\mathrm{e}\Vert_\infty / \Vert\Phi_\mathrm{e}\Vert_\infty$")
    ax1.set_xlabel(r"\# cells per patch per direction")
    ax1.axis(xmin=1e1,xmax=1e4,ymin=3e-7,ymax=1e-1)

    ms = 9
    ax1.loglog(single_ncells,single_Linf_max,'go',ms=ms,mew=0.0,label='single patch')
    ax1.loglog(multi_equiv_ncells3,multi_Linf_max3,'bs',ms=ms,mew=0.0,label='$3^2$ patches')
    ax1.loglog(multi_equiv_ncells9,multi_Linf_max9,'cs',ms=ms,mew=0.0,label='$9^2$ patches')
    ax1.loglog(nested_2l_equiv_ncells,nested_2l_Linf_max,'r^',ms=ms,mew=0.0,label='nested, $l = 2$')
    ax1.loglog(nested_3l_equiv_ncells,nested_3l_Linf_max,'m^',ms=ms,mew=0.0,label='nested, $l = 3$')
    ax1.loglog(xcells,offset*xcells**(-2),'k-',lw=1.0)
    ax1.annotate(xy=(1e2,5.5e-5),s=r"$\propto 1/N^2$",fontsize=16,rotation=-36.5)

    ax1.legend(loc='upper right',frameon=True,numpoints=1,handletextpad=0.0,fontsize=17,fancybox=True)

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
