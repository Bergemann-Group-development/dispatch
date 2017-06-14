# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib import ticker
import dispatch_utils as du
import plot_utils as pu
import numpy as np
import itertools

if __name__ == "__main__":
    
    # set font style
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=15)
    plt.interactive(True)

    offset = 0.0
    data = ['data/input_Ohm']
    run = '.'
    iout = [20]
    npatches = 5
    lsty = ['-','--','-.']
    lcolr = ['b','g','r','c','m','y']
    lcolr = itertools.cycle(lcolr)
    lw = 2.0

    plot_ref_data = False
    if plot_ref_data: ref = np.loadtxt('comp_Masson.txt')

    fig = plt.figure('C-shocks',facecolor='w')
    fig.clf()
    fig.subplots_adjust(left=0.09,right=0.97,bottom=0.10,top=0.97,hspace=0.25,wspace=0.31)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2,sharex=ax1)
    ax3 = fig.add_subplot(2,2,3,sharex=ax1)
    ax4 = fig.add_subplot(2,2,4,sharex=ax1)
    axs = [ax1,ax2,ax3,ax4]
    # Aesthetic tweaks.
    for cax in axs:
        cax.xaxis.set_ticks_position('bottom')
        cax.yaxis.set_ticks_position('left')
        cax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        cax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        for loc,spine in cax.spines.items():
            if loc in ['right','top']:
                spine.set_color('none')

    try:
        iout.index(0)
    except:
        iout.insert(0,0)
    for j in xrange(len(data)):
        for i in iout:
            pp = []
            llcolr = lcolr.next()
            for k in range(npatches):
                p = du.patch(id=k+1,iout=i,run=run,data=data[j])
                pp.append(p)
                li = p.li[0]-1
                ui = p.ui[0]

                if j == 0:
                    ax1.set_ylabel(r'$\rho$ in code units')
                    ax1.axis(ymin=0.36,ymax=1.1,xmin=0.0,xmax=1.0)
                    ax1.plot([p.x[li],p.x[li]],[-100.,100.],'k:',lw=1.0)
                l, = ax1.plot(p.x[li:ui]+offset,p.data[:,0,0,0],linestyle=lsty[j],lw=lw,color=llcolr)
                if lsty[j] == '-.':
                    l.set_dashes([10,3,2,3])
                if plot_ref_data == True:
                    ax1.plot(ref[:,0],ref[:,5])		

                if j == 0:
                    ax2.set_ylabel(r'$B_y$ in code units')
                    ax2.axis(ymin=0.5,ymax=4.0)
                    ax2.plot([p.x[li],p.x[li]],[-100.,100.],'k:',lw=1.0)
                l, = ax2.plot(p.x[li:ui]+offset,p.data[:,0,0,6],linestyle=lsty[j],lw=lw,color=llcolr)
                if lsty[j] == '-.':
                    l.set_dashes([10,3,2,3])
                if plot_ref_data == True:
                    ax2.plot(ref[:,0],ref[:,1])

                if j == 0:
                    ax3.set_xlabel(r'$x$ in code units')
                    ax3.set_ylabel(r'$v_x$ in code units')
                    ax3.axis(ymin=1.5,ymax=5.5)
                    ax3.plot([p.x[li],p.x[li]],[-100.,100.],'k:',lw=1.0)
                l, = ax3.plot(p.x[li:ui]+offset,p.data[:,0,0,2]/p.data[:,0,0,0],linestyle=lsty[j],
                              lw=lw,color=llcolr)
                if lsty[j] == '-.':
                    l.set_dashes([10,3,2,3])
                if plot_ref_data == True:
                    ax3.plot(ref[:,0],ref[:,3])

                if j == 0:
                    ax4.set_xlabel(r'$x$ in code units')
                    ax4.set_ylabel(r'$v_y$ in code units')
                    ax4.axis(ymin=-0.1,ymax=1.6)
                    ax4.plot([p.x[li],p.x[li]],[-100.,100.],'k:',lw=1.0)
                l, = ax4.plot(p.x[li:ui]+offset,p.data[:,0,0,3]/p.data[:,0,0,0],linestyle=lsty[j],
                              lw=lw,color=llcolr)
                if lsty[j] == '-.':
                    l.set_dashes([10,3,2,3])
                if plot_ref_data == True:
                    ax4.plot(ref[:,0],ref[:,4])	

    plt.savefig('C-shock_all.pdf',transparent=True)
    plt.show()
