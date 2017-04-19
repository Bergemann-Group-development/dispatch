import matplotlib.pyplot as pl

def image(f,**kwargs):
    '''
    Turn images the Fortran way; 1st index increases to the right,
    2nd index upwards
    '''
    pl.imshow(f.squeeze().transpose(),origin='lower',**kwargs)
