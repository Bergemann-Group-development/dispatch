m_earth=5.972e27
r_earth=6.371e8
grav=6.673e-8

def scaling():
    class s:
        l=r_earth
        m=m_earth #*1e-10
        v=1e5
        t=l/v
        d=m/l**3
        p=d*v**2
        g=grav*d*t**2
    return s