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
        l=1.
        m=1. #*1e-10
        v=100.
        d=1.
        p=1.
        g=1.
        u=100.
        t=l/u
        e=1.
        temp = 1.
    return s
