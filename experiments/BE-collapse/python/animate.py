import time
import numpy as np
from scipy import special
from numpy import fromfile, float32, prod, size
from matplotlib.pyplot import imshow, show, figure, colorbar
import matplotlib.animation as animation

# Read the data (you may want to have a for loop here)
n=64
shape=(n,n,n)
taskid=1

def read(filename, shape):
    fd=open(filename,'rb')
    dd=fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
    dd=np.transpose(dd)
    im = dd[:,:,n/2]
    im = np.log10(dd[:,:,n/2])
    print 'Density:', np.min(dd), np.max(dd)
    return im

def update(*args):
    global image, i
    i = i + 1
    if(i>27): return
    filename='data/BE_{:0=3d}_{:0=4d}.dat'.format(taskid,i)
    im = read(filename, shape)
    image.set_data(im)

i=0
filename='data/BE_{:0=3d}_{:0=4d}.dat'.format(taskid,i)
fig = figure()
im = read(filename, shape)
# Show a central slice of log(density)
image = imshow(im,interpolation='nearest')
colorbar()
ani = animation.FuncAnimation(fig, update, interval=200, blit=False)
show()
