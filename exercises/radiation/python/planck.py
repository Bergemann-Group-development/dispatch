'''
   Planck function, returns lambda*B_lambda = contribution per logarithmic interval
   Try something similar to this:
       temp=100.                                        # temperature in K
       log_lambda_0=log10(9e-4)                         # lower limit microns
       lamb=logspace(log_lambda_0,log_lamda_0+2.0,9)    # log scale, factor 100
       loglog(lamb,planck(temp,lamb))                   # loglog plot
'''
import numpy as np
from cgs import cgs

def planck(T,lamb):
    h=6.626e-27
    k=1.38064e-16
    c=2.9997e10
    #print h*c/(4e-3*k*T)
    return 2.*h*c**2/lamb**4/(np.exp(h*c/(lamb*k*T))-1.0)