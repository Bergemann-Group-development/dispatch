""" Demo: read and show binary data """
from struct import unpack
from numpy import reshape
from matplotlib.pyplot import imshow, show

s = open("tmp.dat", "rb").read()
data = unpack("@40960f", s)
data = reshape(data,(64,64,10))
for i in range(10):
    print(data[1,1,i])
