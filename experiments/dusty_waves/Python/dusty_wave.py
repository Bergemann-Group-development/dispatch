import numpy as np
#import matplotlib.pyplot as pl


#test = 1
def dusty_wave( x , t , vx_gas , vx_dust , rho_gas , rho_dust , K , cs , k ):
    tg = rho_gas/K
    td = rho_dust/K
    
    a  = -( tg**(-1) + td**(-1) )
    b  =  (k**2)*(cs**2)
    c  = -(k**2)*(cs**2)/td
    
    P  = ((a**2)-3*b)/9.0
    Q  = (a*b/6.0) - (c/2.0) - ((a**3)/27.0)
    D  = (P**3)-(Q**2)
    
    if P<0:
        S  = (1.0/3.0)*arcsinh(Q/np.sqrt(-P**3))
        r0 = (-a/3.0) + 2.0*np.sqrt(-P)*np.sinh(S)
        r1 = (-a/3.0) - np.sqrt(-P)*np.sinh(S) + 1.0j*np.sqrt(-3.0*P)*np.cosh(S)
        r2 = (-a/3.0) - np.sqrt(-P)*np.sinh(S) - 1.0j*np.sqrt(-3.0*P)*np.cosh(S)
    elif P>0 and D<0:
        S  = np.cbrt(Q-np.sqrt(-D))
        R  = np.cbrt(Q+np.sqrt(-D))
        r0 = (-a/3.0) + S + R
        r1 = (-a/3.0) - ((S+R)/2.0) + 1.0j*np.sqrt(3.0)*((S+R)/2.0)
        r2 = (-a/3.0) - ((S+R)/2.0) - 1.0j*np.sqrt(3.0)*((S+R)/2.0)
    elif P>0 and D>0:
        S  = np.arccos(Q/np.sqrt(P**3))
        r0 = (-a/3.0) + 2.0*np.sqrt(P)*np.cos((2.0*np.pi*0.0+S)/3.0)
        r1 = (-a/3.0) + 2.0*np.sqrt(P)*np.cos((2.0*np.pi*1.0+S)/3.0)
        r2 = (-a/3.0) + 2.0*np.sqrt(P)*np.cos((2.0*np.pi*(-1.0)+S)/3.0)
    
    w0 = -1.0*r0
    w1 = -1.0*r1
    w2 = -1.0*r2
    
    # I now need to solve for the remaining 16 coefficients