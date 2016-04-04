from read_data import read_data
from matplotlib.pyplot import loglog,clf,figure,hold,xlabel,ylabel,xlim,ylim
from matplotlib.mlab import find
from scaling import scaling

figure(1)
clf()
hold(1)
n=8
s=scaling()

for i in range(2,n+1):
    f=read_data(15+27*i,0)
    loglog(f.x[:],f.d[:,12,12]*s.d,'-+')

f=read_data(15+27*n-1,0)

w=find(f.x>1)
w=w[0]
f.d[13:w,12,12]=10.*f.d[w,12,12]
#loglog(f.x[13:],f.d[13:,12,12]*s.d,'-+')
ylim(1e-12,2e-8)
xlim(0.5,5e3)
ylabel('density [cgs]')
xlabel('x [Earth radii]')