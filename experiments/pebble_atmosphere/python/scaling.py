m_earth=5.972e27
r_earth=6.371e8
grav=6.673e-8
yr = 3.156e+7
au = 1.498e+13
kms = 1e+5
mu = 2.
m_u = 1.6726e-24
k_b = 1.3807e-16

def scaling():
    class s:
        l=r_earth
        m=m_earth #*1e-10
        v=1e5
        t=l/v
        d=m/l**3
        p=d*v**2
        g=grav*d*t**2
        u=kms
        e=m*u**2
        temp = mu*(m_u/m)/(k_b/e)
    return s
