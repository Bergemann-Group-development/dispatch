from read_data import read_data
from matplotlib.pyplot import loglog,clf,figure,hold,xlabel,ylabel,ylim
from matplotlib.mlab import find
from scaling import scaling

figure(1)
clf()
hold(1)
n=8
s=scaling()

for i in range(2,n+1):
    f=read_data(23+27*i,0)
    loglog(f.z[:],f.d[12,12,:]*s.d,'-+')

f=read_data(23+27*n-9,0)

w=find(f.z>1)
w=w[0]-1
f.d[12,12,13:w]=f.d[12,12,w]
loglog(f.z[13:],f.d[12,12,13:]*s.d,'-+')
ylim(1e-14,3e-8)
ylabel('density [cgs]')
xlabel('z [Earth radii]')