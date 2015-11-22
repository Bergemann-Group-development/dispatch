""" read and show advection data """
from numpy import reshape, fromfile, float32, prod
from matplotlib.pyplot import imshow, show

shape = (64,64)
fd = open('snapshot.dat', 'rb')

while True:
  data = fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)
#print(data)
  imshow(data,origin='lower')
  show()
