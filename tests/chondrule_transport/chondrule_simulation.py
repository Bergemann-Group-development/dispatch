import numpy as np
import matplotlib.pyplot as pl

# Chondrule object, with attributes
class chondrule:
    age=0.0
    nreset=0
    lost=False

# Make a list with n of them    
n=300000
pp=[]
for i in range(n):
    pp.append(chondrule())

# Setup time simulation    
time=0.0
dtime=0.05          # time step
p_loss=2.1          # probability of loss per Myr
p_reset=0.30        # probability of reset per Myr
p_accrete=0.0001    # probability of accretion per Myr

# Loop over time, accumulate accreted chondrules
accr=[]
while time<4.2:
    time=time+dtime
    for p in pp:
        ran=np.random.random(3)
        # With probability p_loss, loose it
        if ran[0] < p_loss*dtime:
            p.lost=True
        # With probability p_reset, reset it
        if ran[1] < p_reset*dtime:
            p.nreset+=1
            p.age=time
        # With probability p_accrete, save it as accreted and loose it
        if ran[2] < p_accrete*dtime:
            p.lost=True
            accr.append(p)
    # collect remaining chondrules (not lost or accreted)
    pp=[p for p in pp if not p.lost]
    # print time, remaining chondrules, accreted chondrules
    print '{:5.2f} {:6d} {:5d}'.format(time,np.size(pp),np.size(accr))

# Accumulate arrays with ages (accreted and remaining)
ages1=[]
for p in accr:
    ages1.append(p.age)
ages2=[]
for p in pp:
    ages2.append(p.age)

# Plot histograms
pl.clf()
pl.figure(1,figsize=(9,7))
pl.subplot(2,1,1)
pl.hist(ages1,bins=5)
pl.title('continuous accretion')
pl.xlabel('Myr')
pl.subplot(2,1,2)
pl.hist(ages2,bins=5)
pl.title('late accretion')
pl.xlabel('Myr')
pl.tight_layout()
pl.show()