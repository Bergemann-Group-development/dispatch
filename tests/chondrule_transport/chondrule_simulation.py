import numpy as np
import matplotlib.pyplot as pl

# Chondrule object, with attributes
class chondrule():
    created=0.     # time created
    age=0.         # time of last reset
    nreset=0       # number of resets
    lost=False     # lost or not

# Rate object, with attributes
class rates():
    case=''        # rate case 
    loss=0.        # rate of loss per Myr
    reset=0.       # rate of reset per Myr
    accrete=0.     # rate of accretion per Myr
    ctime=0.       # rate of production decay per Myr

# Initial number of chondrules
n=1000000/0.3

# Setup time simulation
dtime=0.0501        # time step
tend=4.2

troels=rates()
troels.case='EXTENDED CREATION'
troels.loss=7.0
troels.reset=25.0
troels.accrete=0.005
troels.ctime=0.5

aake=rates()
aake.case='BRIEF CREATION'
aake.loss=2.0
aake.reset=0.25
aake.accrete=0.0005
aake.ctime=0.2

# Loop over cases: Figures can be modified interactive, before (re-)saving as .png
cases=[aake,troels]
for r in cases:
    # Loop over time, accumulate accreted chondrules
    pp=[]
    accr=[]
    nold=0.
    time=0.1    
    while time<tend:
        time=time+dtime
        for p in pp:
            ran=np.random.random(3)
            # With probability r.loss, loose it
            if ran[0] < r.loss*dtime:
                p.lost=True
            # With probability r.reset, reset it
            pp_reset = r.reset*dtime
            #if (time > 2.5):
            #pp_reset = pp_reset * np.exp(-20.*(time-2.5))
            if ran[1] < pp_reset:
                p.nreset+=1
                p.age=time
            # With probability r.accrete, save it as accreted and loose it
            if ran[2] < r.accrete*dtime:
                p.lost=True
                accr.append(p)
        # collect remaining chondrules (not lost or accreted)
        pp=[p for p in pp if not p.lost]
    
        # add new chondrules, with exponentially decaying rate
        dndt = n * np.exp(-time / r.ctime)
        nnew = nold + dndt * dtime
        nstep = int(np.floor(nnew) - np.floor(nold))
        nold = nnew
        for i in range(nstep):
            pp.append(chondrule())
            pp[-1].age = time
            pp[-1].create = time
    
        # print time, remaining chondrules, accreted chondrules, new chondrules
        print 't={:5.2f}  ncur={:6d}  nacc={:5d}  nnew={:5d}'.format(time,np.size(pp),np.size(accr),nstep)
    
    # Accumulate arrays with ages (accreted and remaining for Troels, remaining for Aake)
    ages=[]
    nreset=[]
    created=[]
    if r==troels:
        all=np.append(pp,accr)
    else:
        all=pp
    for p in all:
        created.append(p.create)
        ages.append(p.age)
        nreset.append(p.nreset)
    
    # Figure size
    pl.figure(figsize=(9,7))
    pl.clf()
    pl.interactive(True)

    bins = np.linspace(0.,tend+dtime,5+1)
    print bins
    
    # age distribution
    pl.subplot(2,2,1)
    pl.hist(ages,bins=bins)
    pl.title('age distribution')
    pl.xlabel('Myr')
    pl.xlim(0.0,tend+dtime)
    
    # age - resetting scatter plot
    pl.subplot(2,2,3)
    pl.scatter(ages, nreset)
    pl.xlabel('Myr')
    pl.ylabel('nreset')
    pl.xlim(0.0,tend+dtime)
    
    # creation time histogrm
    pl.subplot(2,2,2)
    pl.hist(created,bins=bins,normed=True)
    pl.title('creation time distribution')
    pl.xlabel('Myr')
    pl.xlim(0.0,tend+dtime)
    
    # evolution time as a function of age
    pl.subplot(2,2,4)
    pl.scatter(ages,np.array(ages)-np.array(created))
    pl.xlabel('Myr')
    pl.ylabel('evolved [Myr]')
    pl.xlim(0.0,tend+dtime)
    
    # Add main title, and save as .png (if necessary, modify margins interactively)
    pl.suptitle(r.case)
    pl.tight_layout()
    pl.show()
    pl.pause(0.001)
    pl.savefig(r.case+'.png')
