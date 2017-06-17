import numpy as np
import matplotlib.pyplot as pl

# Chondrule object, with attributes
class chondrule():
    created=0.
    age=0.
    nreset=0
    lost=False
# Rate attributes
class rates():
    case=''        # rate case 
    loss=0.        # rate of loss per Myr
    reset=0.       # rate of reset per Myr
    accrete=0.     # rate of accretion per Myr
    ctime=0.       # rate of production decay per Myr

# Chondrules produced per Myr
# Assumed to decay exponentially
def create_chondrule(time, n, ctime):
    nc = n * np.exp(-time / ctime)
    return nc

# Initial number of chondrules
n=1000000/0.3

# Setup time simulation
dtime=0.0201        # time step
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
            # With probability p_loss, loose it
            if ran[0] < r.loss*dtime:
                p.lost=True
            # With probability p_reset, reset it
            pp_reset = r.reset*dtime
            #if (time > 2.5):
            #pp_reset = pp_reset * np.exp(-20.*(time-2.5))
            if ran[1] < pp_reset:
                p.nreset+=1
                p.age=time
                #if (time > 2.5):
                #    print 'Time reset for chondrule :', time, p.nreset
            # With probability p_accrete, save it as accreted and loose it
            if ran[2] < r.accrete*dtime:
                p.lost=True
                accr.append(p)
        # collect remaining chondrules (not lost or accreted)
        pp=[p for p in pp if not p.lost]
    
        # add new chondrules
        dndt = create_chondrule(time,n,r.ctime)
        nnew = nold + dndt * dtime
        nstep = int(np.floor(nnew) - np.floor(nold))
        nold = nnew
        for i in range(nstep):
            pp.append(chondrule())
            pp[-1].age = time
            pp[-1].create = time
    
        # print time, remaining chondrules, accreted chondrules
        print 't={:5.2f} ncur={:6d} nacc={:5d} nnew={:5d}'.format(time,np.size(pp),np.size(accr),nstep)
    
    # Accumulate arrays with ages (accreted and remaining)
    ages1=[]
    nreset=[]
    created=[]
    if r==troels:
        all=np.append(pp,accr)
    if r==aake:
        all=pp
    
    for p in all:
        created.append(p.create)
        ages1.append(p.age)
        nreset.append(p.nreset)
    
    # Plot histograms
    pl.figure(figsize=(9,7))
    pl.clf()
        
    bins = np.linspace(0.,tend+dtime,5+1)
    print bins
    
    pl.subplot(2,2,1)
    pl.hist(ages1,bins=bins)
    pl.title('age distribution')
    pl.xlabel('Myr')
    pl.xlim(0.0,tend+dtime)
    
    pl.subplot(2,2,3)
    pl.scatter(ages1, nreset)
    pl.xlabel('Myr')
    pl.ylabel('nreset')
    pl.xlim(0.0,tend+dtime)
    
    pl.subplot(2,2,2)
    pl.hist(created,bins=bins,normed=True)
    pl.title('creation time distribution')
    pl.xlabel('Myr')
    pl.xlim(0.0,tend+dtime)
    
    pl.subplot(2,2,4)
    pl.scatter(ages1,np.array(ages1)-np.array(created))
    pl.xlabel('Myr')
    pl.ylabel('evolved [Myr]')
    pl.xlim(0.0,tend+dtime)
    
    pl.suptitle(r.case)
    pl.tight_layout()
    pl.draw()
    pl.show()
    pl.savefig(r.case+'.png')
