""" read and show advection data """
from numpy import reshape, fromfile, float32, prod, concatenate, log
from matplotlib.pyplot import imshow, show, figure, subplot

shape = (100,100,100)

fd = open('objects_01.dat', 'rb')
data1 = fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)

fd = open('objects_02.dat', 'rb')
data2 = fromfile(file=fd, dtype=float32, count=prod(shape)).reshape(shape)

fig = figure()
subplot(211)
imshow(data1[1:200,1:100,50],origin='lower')
subplot(212)
imshow(data2[1:200,1:100,50],origin='lower')
show()
