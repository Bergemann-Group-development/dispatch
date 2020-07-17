import scipy
import numpy as np

def interp(x,y,x0):
    """ Interpolate in x[:] and y[:] to x-coordinat(s) x0
    """
    kind='cubic'
    f=scipy.interpolate.interp1d(x,y,kind,fill_value='extrapolate')
    return f(x0)

def deriv(x,y):
    return np.gradient(y,x)

def integral(x,y,y1=None):
    """ Integrate over x[:],y[:] with optional derivate y1[:] 
    """
    if y1 is None:
        y1=deriv(x,y)
    f=np.zeros(len(x))
    for i in range(1,len(x)):
        dx=x[i]-x[i-1]
        f[i]=f[i-1]+dx*(y[i]+y[i-1])/2+dx**2*(y1[i-1]-y1[i])/12
    return f