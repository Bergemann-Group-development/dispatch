import scipy
import numpy as np

def interp(x,y,x0):
    """ Interpolate in x[:] and y[:] to x-coordinat(s) x0
    """
    kind='cubic'
    f=scipy.interpolate.interp1d(x,y,kind,fill_value='extrapolate')
    return f(x0)

def deriv(x,y=None,edge_order=2,**kw):
    if y is None:
        return np.gradient(x,edge_order=edge_order,**kw)
    else:
        return np.gradient(y,x,edge_order=edge_order,**kw)

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

def open(case='stagger',file=None,top='~/codes/dispatch2',verbose=0):
    import dispatch.EOS
    import os
    if case=='stagger':
        os.chdir(top+'/data/eos/stagger')
        eos=dispatch.EOS.stagger
    elif case=='tomida':
        os.chdir(top+'/data/eos/tomida')
        eos=dispatch.EOS.tomida
    elif case=='bifrost':
        os.chdir(top+'/data/eos/bifrost')
        eos=dispatch.EOS.bifrost
    return eos(verbose=verbose) if file is None else eos(file,verbose=verbose)
