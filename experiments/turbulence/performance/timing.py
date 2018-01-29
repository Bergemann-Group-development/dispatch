import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
import itertools

if __name__ == "__main__":
    #plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('font', size=17)
    plt.interactive(True)
    fig1 = plt.figure('timing',facecolor='w',figsize=(5.4,3.3))
    
    ms = 5.0

    # timing data
    # RAMSES
    #ramses_nnodes = np.array([1, 8, 64, 512])
    #ramses_ncores = np.array([      1,      3,      6,     12,     24])
    #ramses_mus    = np.array([[ 19.07,  19.67,  20.74,  27.54,  41.48],
                              #[ 20.5 ,  21.46,  25.75,  33.9 ,  49.35],
                              #[ 20.5 ,  21.82,  26.46,  35.48,  48.78],
                              #[  0.  ,   0.  ,  21.74,  36.62,  48.07]])

    # DISPATCH/HLRS/Hazel Hen/HLLC
    hlrs_ranks  = np.array([1, 8, 64, 256, 512, 1024, 1000, 2000])
    hlrs_ncores = np.array( [  1   ,  6   , 12   , 24   ])
    hlrs_mus    = np.array([[  0.29,  0.43,  0.60,  0.62],
                            [  0.38,  0.54,  0.64,  0.63],
                            [  0.39,  0.50,  0.63,  0.68],
                            [  0.56,  0.57,  0.65,  0.69],
                            [  0.57,  0.57,  0.70,  0.71],
                            [  0.57,  0.59,  0.72,np.nan],
                            [np.nan,np.nan,np.nan,np.nan],
                            [np.nan,np.nan,np.nan,np.nan] ])

    c = 24
    hlrs_cores  = np.array([1, 6, 12, 24, 2*c, 4*c, 8*c, 16*c, 32*c, 64*c, 128*c, 256*c, 512*c, 1024*c, 1000*c, 2000*c])
    hlrs_ncells = np.array( [   64])
    hlrs_mus2   = np.array([[  0.29], # 1 core
                            [  0.43], # 6 cores
                            [  0.60], # 12 cores
                            [  0.62], # 24 cores
                            [  0.65], # 2 nodes (48 cores)
                            [  0.64], # 4 nodes (96 cores)
                            [  0.63], # 8 nodes (192 cores)
                            [  0.63], # 16 nodes (384 cores)
                            [  0.63], # 32 nodes (768 cores)
                            [  0.68], # 64 nodes (1536 cores)
                            [  0.65], # 128 nodes (3072 cores)
                            [  0.69], # 256 nodes (6144 cores)
                            [  0.71], # 512 nodes (12288 cores)
                            [  0.76], # 1024 nodes (24576 cores)
                            [np.nan], # 1000 nodes
                            [np.nan]]) # 2000 nodes

    # DISPATCH/Marconi/KNL/HLLC
    mknl_ranks  = np.array([1, 8, 64, 256, 512, 1024, 1000, 2000])
    mknl_ncores = np.array( [  1,  2,     4,     8,    17,    34,    68])
    # using 64**3 cells/patch and 8**3 patches/rank
    mknl_mus    = np.array([[  0.96,  0.97,  0.96,  0.96,  0.96,  0.97,  0.98], # 1 node
                            [  0.97,np.nan,  0.99,  0.99,  0.99,  1.07,  1.09], # 8 nodes
                            [  1.00,np.nan,  1.01,  0.99,  1.00,  1.09,  1.11], # 64 nodes
                            [  0.99,np.nan,  0.99,  0.98,  0.98,  1.07,  1.11], # 256 nodes
                            [  0.98,  1.03,  0.98,  0.98,  0.99,  1.07,  1.09], # 512 nodes
                            [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan], # 1024 nodes
                            [  0.98,  1.01,  0.99,  0.98,  1.00,  1.05,  1.11], # 1000 nodes
                            [  1.01,np.nan,  1.01,  0.98,  0.99,  1.08,np.nan] ]) # 2000 nodes

    # using 32*83 cells/patch and 8**3 patches/rank
    mknl_mus32  = np.array( [  1.31,  1.28,  1.23,  1.31,  1.44,  1.63,  1.60])

    # now with varying number of cells per patch
    mknl_ncells = np.array( [    32,    48,    56,    64,    72])
    mknl_mus2   = np.array([[  1.35,  1.02,  0.98,  0.96,  0.92], # 1 core
                            [  1.31,  1.02,  0.98,  0.97,  0.94], # 2 cores
                            [  1.21,  1.03,  0.99,  0.96,  0.92], # 4 cores
                            [  1.20,  1.01,  0.98,  0.96,  0.91], # 8 cores
                            [  1.25,  1.02,  0.99,  0.96,  0.91], # 17 cores
                            [  1.28,  1.06,  1.01,  0.97,  0.93], # 34 cores
                            [  1.34,  1.04,  1.03,  0.98,  0.95], # 68 cores
                            [  1.34,  1.10,  1.05,  1.03,  0.99], # 4 ranks
                            [  1.55,  1.20,  1.12,  1.09,  1.20], # 8 ranks
                            [  1.58,  1.20,  1.14,  1.09,np.nan], # 16 ranks
                            [  1.56,  1.18,  1.12,  1.08,np.nan], # 32 ranks
                            [  1.63,  1.24,  1.14,  1.11,  1.27], # 64 ranks
                            [  1.59,  1.22,  1.14,  1.09,np.nan], # 128 ranks
                            [  1.59,  1.20,  1.14,  1.11,  1.15], # 256 ranks
                            [  1.60,  1.27,  1.15,  1.09,  1.29], # 512 ranks
                            [np.nan,np.nan,np.nan,np.nan,np.nan], # 1024 ranks
                            [  1.60,  1.19,  1.13,  1.11,  1.27], # 1000 ranks
                            [  1.60,  1.20,  1.13,  1.08,np.nan] ]) # 2000 ranks
    c = 34
    mknl_cores = np.array([1, 2, 4, 8, 17, 34, 68, 4*c, 8*c, 16*c, 32*c, 64*c,
                           128*c, 256*c, 512*c, 1024*c, 1000*c, 2000*c])

    # DISPATCH/Marconi/KNL/Staggger2e
    mknl_ncells_stagger = np.array( [    32,    48,    56,    64,    72])
    mknl_mus2_stagger   = np.array([[], # 1 core
                                    [], # 2 cores
                                    [], # 4 cores
                                    [], # 8 cores
                                    [], # 17 cores
                                    [], # 34 cores
                                    [], # 68 cores
                                    [], # 8 ranks
                                    [], # 64 ranks
                                    [] ]) # 256 ranks
    c = 34
    mknl_cores_stagger = np.array([1, 2, 4, 8, 17, 34, 68, 8*c, 64*c, 256*c])

    hlrs_ctot = hlrs_ranks[:,np.newaxis]*hlrs_ncores[np.newaxis,:]
    mknl_ctot = mknl_ranks[:,np.newaxis]*mknl_ncores[np.newaxis,:]

    mknl_ncells_str = [r'{0}$^3$ HLLC'.format(m) for m in mknl_ncells]

    # End data; Start plotting.
    lcolrs = ['r','b','g','c','m','Orange','y','Navy']
    lcolrs1 = itertools.cycle(lcolrs)
    lcolrs2 = itertools.cycle(lcolrs)
    lcolrs3 = itertools.cycle(lcolrs)

    fig1.clf()
    fig1.subplots_adjust(left=0.14,bottom=0.16,right=0.87,top=0.97)
    ax1 = fig1.add_subplot(1,1,1)

    # RAMSES
    #for i in range(ramses_mus.shape[0]):
    #    ax1.semilogy(ramses_ncores[0:6],ramses_mus[i,0:6],'o-',color=lcolrs1.next(),
    #               mew=0.0,lw=2.0)
    #plt.figtext(0.19,0.69,'RAMSES unigrid, 1-512 ranks',size=sz)

    # DISPATCH/HazelHen/Haswell
    #ax1.loglog(hlrs_ctot[0,:],hlrs_mus[0,:],'o-',color=lcolrs2.next(),
               #mew=0.0,label=hlrs_ranks[i],lw=2.0,ms=ms)
    #for i in range(1,hlrs_mus.shape[0]):
        #ax1.loglog(hlrs_ctot[i,3],hlrs_mus[i,3],'o-',color=lcolrs2.next(),
                   #mew=0.0,label=hlrs_ranks[i],lw=2.0,ms=ms)
    #plt.figtext(0.47,0.34,'Xeon Haswell',size=11.0)

    # DISPATCH/Marconi/KNL
    #for i in range(mknl_mus.shape[0]):
        #l, = ax1.loglog(mknl_ctot[i,:],mknl_mus[i,:],'o-',color=lcolrs3.next(),
                        #mew=0.0,lw=2.0,ms=ms)
        ##l.set_dashes([2,2])
    #plt.figtext(0.34,0.55,'Xeon Phi KNL',size=11.0)

    # DISPATCH/HazelHen/Haswell
    for i in range(hlrs_mus2.shape[1]):
        ax1.loglog(hlrs_cores,hlrs_mus2[:,i],'s',color='g',
                   mew=0.0,lw=1.0,ms=ms)
    plt.figtext(0.29,0.40,'Xeon Haswell',size=11.0)

    # DISPATCH/Marconi/Knights Landing
    for i in [0,1,3]:
        ax1.loglog(mknl_cores,mknl_mus2[:,i],'o',color=lcolrs3.next(),
                   mew=0.0,lw=1.0,ms=ms,label=mknl_ncells_str[i])

    # DISPATCH/Marconi/KNL/Stagger
    n=13
    y32=np.ones(n)
    y40=np.ones(n)
    y48=np.ones(n)
    c  =np.zeros(n)
    i=0  ; y32[i]=1.05; y40[i]=.88; y48[i]=.75; c[i]=4      # 4 core
    i=i+1; y32[i]=1.07; y40[i]=.90; y48[i]=.75; c[i]=8      # 8 cores
    i=i+1; y32[i]=1.15; y40[i]=.85; y48[i]=.75; c[i]=17     # 17 cores
    i=i+1; y32[i]=1.75; y40[i]=.85; y48[i]=.75; c[i]=34     # 34 cores
    i=i+1; y32[i]=0.97; y40[i]=.85; y48[i]=.75; c[i]=68     # 1 node
    i=i+1; y32[i]=0.99; y40[i]=.82; y48[i]=.77; c[i]=68*2   # 2 nodes
    i=i+1; y32[i]=0.99; y40[i]=.83; y48[i]=.80; c[i]=68*4   # 4 nodes
    i=i+1; y32[i]=1.00; y40[i]=.83; y48[i]=.80; c[i]=68*8   # 8 nodes
    i=i+1; y32[i]=1.01; y40[i]=.87; y48[i]=.81; c[i]=68*16  # 16 nodes
    i=i+1; y32[i]=1.05; y40[i]=.85; y48[i]=.82; c[i]=68*32  # 32 nodes
    i=i+1; y32[i]=1.03; y40[i]=.85; y48[i]=.80; c[i]=68*64  # 64 nodes
    i=i+1; y32[i]=0.95; y40[i]=.82; y48[i]=.80; c[i]=68*128 # 128 nodes
    i=i+1; y32[i]=0.90; y40[i]=.81; y48[i]=.79; c[i]=68*256 # 256 nodes
                   
    # DISPATCH/Marconi/KNL/Stagger
    ax1.loglog(c,y48,'^',color=lcolrs3.next(),
                mew=1.0,lw=1.0,ms=ms,label=r'48$^3$ STG')
    plt.figtext(0.22,0.64,'Xeon Phi KNL',size=11.0)

    ax1.set_xlabel(r'cores',fontsize=15)
    ax1.set_ylabel(r'core-$\mu$s/cell-update',fontsize=15,va='center')

    ax1.tick_params(axis='both', which='major', labelsize=12)

    leg = ax1.legend(loc='upper left',frameon=False,numpoints=1,fontsize=11,
                     borderaxespad=0.35,handletextpad=0.3,handlelength=1.5,ncol=2,
                     columnspacing=0.60)
    #leg.set_title('MPI ranks',prop={'size':11})

    ax1.axis(ymin=0.1,ymax=10.0)

    # parasite axis for update time
    ax2 = ax1.twinx()
    ax2.set_ylabel('$2000^3$ update time (s)',fontsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=12)

    # Enter strong scaling data here:
    m_strong=128
    ppt=1e-6*np.array([1.1,0.9,0.85,0.90,1.1])
    x=68*2**np.arange(4,9)
    time=ppt*m_strong**3*x[2]/x
    
    # Perfect scaling line
    x1=68*2**np.arange(3.5,9.2)
    time1=ppt.min()*m_strong**3*x[2]/x1
    ax2.loglog(x1,time1)
    ax2.loglog(x,time,'^')

    ax2.axis(xmin=0.7,xmax=2.0e5)
    ax2.axis(ymin=1.0/8,ymax=2.21e2/8)

    # Aesthetic tweaks.
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    #ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    for loc,spine in ax1.spines.items():
        if loc in ['top','right']:
            spine.set_color('none')
    for loc,spine in ax2.spines.items():
        if loc in ['top']:
            spine.set_color('none')

    plt.draw()
    plt.savefig('timing.pdf',transparent=True)
