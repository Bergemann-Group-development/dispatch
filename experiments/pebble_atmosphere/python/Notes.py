import mayavi.mlab as ml
import numpy as np
'''
    show_traces('pebbles6',radius=50)
        shows fast particles entering the Hill sphere
    show_traces('pebbles8',radius=200)
        shows slow particles, both initial and drifting in
'''

def cube(x,y,z,sz,phi):
    c=np.cos(phi+np.pi/4)
    s=np.sin(phi+np.pi/4)
    x1=0.5*sz*c*np.sqrt(2.)
    y1=0.5*sz*s*np.sqrt(2.)
    z1=0.5*sz
    x2=x+np.array([+x1,-y1,-x1,+y1,+x1])
    y2=y+np.array([+y1,+x1,-y1,-x1,+y1])
    z2=z+np.array([+z1,+z1,+z1,+z1,+z1])
    lw=0.1
    ml.plot3d(x2,y2,+z2,tube_radius=None,line_width=lw)
    ml.plot3d(x2,y2,-z2,tube_radius=None,line_width=lw)
    for i in range(4):
        ml.plot3d([x2[i],x2[i]],[y2[i],y2[i]],[-z1,z1],tube_radius=None,line_width=lw)

au=1.5e13
re=6.371e8
szp=2.2*re*0.9**(1./3.)
n_level=7
szc=szp*3.**n_level
f=2.*np.pi*au/szc
print f

ml.clf()

nphi=31
radius=szc*31/(2.*np.pi)
dphi=2.*np.pi/nphi
for i in range(nphi):
    phi=i*dphi+np.pi/4.
    x=radius*np.cos(phi)
    y=radius*(np.sin(phi)-1.0)
    z=0.0
    cube(x,y,z,szc,phi)
    #print i,x/szc
    
ml.view(focalpoint=[0.,0.,0.],azimuth=90.,elevation=60.)