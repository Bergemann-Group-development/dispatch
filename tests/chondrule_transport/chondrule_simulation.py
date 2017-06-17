import numpy as np
import matplotlib.pyplot as pl

# Chondrule object, with attributes
class chondrule():
    age=0.
    nreset=0
    lost=False

# Chondrules produced per Myr
# Assumed to decay exponentially
def create_chondrule(time, n):
    nc = n * np.exp(-time / ctime)
    return nc

# Initial number of chondrules
n=1000000/0.3

# Make a list with n of them
pp=[]
#for i in range(n):
#    pp.append(chondrule())

pp.append(chondrule())

# Setup time simulation
time=0.1
dtime=0.0101        # time step
tend=4.2

p_loss=7.0         # probability of loss per Myr
p_reset=25.0        # probability of reset per Myr
p_accrete=0.005  # probability of accretion per Myr
ctime = 0.5         # time decay constant for fresh chondrules

# Loop over time, accumulate accreted chondrules
accr=[]
nold=0.

while time<tend:
    time=time+dtime
    for p in pp:
        ran=np.random.random(3)
        # With probability p_loss, loose it
        if ran[0] < p_loss*dtime:
            p.lost=True
        # With probability p_reset, reset it
        pp_reset = p_reset*dtime
        #if (time > 2.5):
          #pp_reset = pp_reset * np.exp(-20.*(time-2.5))
        if ran[1] < pp_reset:
            p.nreset+=1
            p.age=time
            #if (time > 2.5):
            #    print 'Time reset for chondrule :', time, p.nreset
        # With probability p_accrete, save it as accreted and loose it
        if ran[2] < p_accrete*dtime:
            p.lost=True
            accr.append(p)
    # collect remaining chondrules (not lost or accreted)
    pp=[p for p in pp if not p.lost]

    # add new chondrules
    dndt = create_chondrule(time,n)
    nnew = nold + dndt * dtime
    nstep = int(np.floor(nnew) - np.floor(nold))
    nold = nnew
    for i in range(nstep):
        pp.append(chondrule())
        pp[-1].age = time

    # print time, remaining chondrules, accreted chondrules
    print 't={:5.2f} ncur={:6d} nacc={:5d} nnew={:5d}'.format(time,np.size(pp),np.size(accr),nstep)

# Accumulate arrays with ages (accreted and remaining)
ages1=[]
nreset=[]
for p in pp:
    ages1.append(p.age)
    nreset.append(p.nreset)

for p in accr:
    ages1.append(p.age)
    nreset.append(p.nreset)

# Plot histograms
pl.clf()

pl.subplot(2,1,1)
bins = np.linspace(0.,tend+dtime,5+1)
print bins
pl.hist(ages1,bins=bins,normed=True)
pl.title('chondrule distribution')
pl.xlabel('Myr')
pl.xlim(0.0,tend+dtime)

pl.subplot(2,1,2)
pl.scatter(ages1, nreset)
#pl.hist(ages2,bins=np.arange(0.,tend,tend/5))
#pl.title('late accretion')
pl.xlabel('Myr')
pl.ylabel('nreset')
pl.xlim(0.0,tend+dtime)

pl.show()
