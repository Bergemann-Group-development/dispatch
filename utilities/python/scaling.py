from __future__ import print_function

class cgs:
    name='cgs'
    m_earth=5.972e27
    m_sun=1.989e33
    r_earth=6.371e8
    grav=6.673e-8
    yr=3.156e+7
    au=1.498e+13
    kms=1e+5
    mu=2.
    m_u=1.6726e-24
    k_b=1.3807e-16
    h_p=6.6260e-27
    e=4.8032e-10
    c=2.9979e10
    stefan=5.6704e-8

class SI:
    name='SI'
    m_earth=5.972e24
    m_sun=1.989e30
    r_earth=6.371e6
    grav=6.673e-11
    yr=3.156e+7
    au=1.498e+11
    kms=1e+3
    mu=2.
    m_u=1.6726e-27
    k_b=1.3807e-23
    h_p=6.62606e-34
    e=1.6022e-19
    c=2.9979e8
    stefan=5.6704e-5

def scaling(type='ISM',units=cgs,verbose=0,mu=2):
    '''
    Return a structure with scaling constants.  Use e.g.

        scgs=scaling(cgs)
        sSI=scaling(SI)
        print (cgs.k_b, SI.k_b)
    '''
    if verbose>0:
        print("using "+units.name+" units")
    class s:
        system=units.name
        if type=='ISM':
            l=units.pc
            d=1e-24
            v=1e5
            t=l/v
        elif type=='solar':
            l=1e8
            t=1e2
            d=1e-7
            v=l/t
        m=d*l**3
        p=d*v**2
        g=units.grav*d*t**2
        u=units.kms
        e=m*u**2
        temp = mu*(units.m_u)/(units.k_b)*v**2
    return s
