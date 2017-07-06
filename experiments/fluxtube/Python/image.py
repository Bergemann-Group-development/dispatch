# -*- coding: utf-8 -*-
'''
   To adhere to normal imaging convention, where the first index
   is horizontal, and the second index increases upwards, and to
   get rid of singleton dimensions in the input array, we need a
   squence of boring operators, which we are better off putting
   away in a file, such as this one
'''
import matplotlib.pyplot as pl
import matplotlib.cm as cm
def image(f,title=''):
    pl.clf()
    pl.title(title)
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest') #, vmin=6.0, vmax=10.0)
    pl.colorbar()
    pl.show(block=False)
