import matplotlib.pyplot as pl
import dispatch_utils as du ; reload(du)
import plot_utils as pu     ; reload(pu)

run='m05r27t10b'

p=du.patch(204,10,run)
p.stripit()

pl.figure(1)
pl.clf()
pl.show()
fmax=0.0
fmin=0.0
for i in range(20):
    pl.subplot(4,5,i+1)
    f=p.u1[:,:,i]
    pu.image(f)
    if f.min() < fmin: fmin=f.min()
    if f.max() > fmax: fmax=f.max()
pl.tight_layout()
print fmin,fmax

pl.figure(2)
pl.clf()
pl.show()
fmax=0.0
fmin=0.0
for i in range(20):
    pl.subplot(4,5,i+1)
    f=p.u2[:,:,i]
    pu.image(f)
    if f.min() < fmin: fmin=f.min()
    if f.max() > fmax: fmax=f.max()

pl.figure(3)
pl.clf()
pl.show()
fmax=0.0
fmin=0.0
for i in range(20):
    pl.subplot(4,5,i+1)
    f=p.u3[:,i,:]
    pu.image(f)
    if f.min() < fmin: fmin=f.min()
    if f.max() > fmax: fmax=f.max()
pl.tight_layout()
print fmin,fmax

pl.figure(4)
pl.clf()
pl.show()
p0=du.patch(204,0,run)
p0.stripit()
fmax=0.0
fmin=0.0
for i in range(20):
    pl.subplot(4,5,i+1)
    f=1.-p.d[:,i,:]/p0.d[:,i,:]
    pu.image(f)
    if f.min() < fmin: fmin=f.min()
    if f.max() > fmax: fmax=f.max()
pl.tight_layout()
print fmin,fmax
