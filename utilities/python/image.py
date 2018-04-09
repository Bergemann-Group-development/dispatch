'''
   To adhere to normal imaging convention, where the first index
   is horizontal, and the second index increases upwards, and to
   get rid of singleton dimensions in the input array, we need a
   squence of boring operators, which we are better off putting
   away in a file, such as this one
'''
import matplotlib.pyplot as pl
def image(f,title='',vlim=[None,None]):
    pl.clf()
    pl.title(title)
    pl.imshow(f.squeeze().transpose(), 
      origin='lower',interpolation='nearest',vmin=vlim[0],vmax=vlim[1])
    pl.colorbar()
    pl.show(block=False)
