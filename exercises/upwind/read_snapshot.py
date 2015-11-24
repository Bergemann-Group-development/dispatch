""" read and show advection data """
from numpy import reshape, fromfile, float32, prod, concatenate, log
from matplotlib.pyplot import imshow, show

shape = (100,100,100)

fd = open('objects_01.dat', 'rb')
data1 = fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)

fd = open('objects_02.dat', 'rb')
data2 = fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)

data = concatenate((data1,data2), axis=0)

imshow(log(data[1:200,1:100,50]),origin='lower')
show()
