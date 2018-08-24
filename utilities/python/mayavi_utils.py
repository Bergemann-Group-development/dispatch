# Pythn 2/3 compatibility
from __future__ import print_function

import matplotlib.mlab as mml
import mayavi.mlab as ml
import numpy as np
import os

def mayavi_utils():
    '''
    spin():       Spin the view interactively
    zoom():       Zoom the view a factor exp(e), in n steps
    dolly():      Dolly the view a factor exp(e), in n steps
    spin_zoom():  Spin + zoom the view interactively
    spin_movie(): Spin one turn, saving images in dir
    cube():       Construct a cube
    torus():      Construct a torus of cubes
    
    Type e.g. spin?<RETURN> for iPython help, or spin( for 
    pop-up help in Canopy.
    '''
    print("Type mayavi_utils?<RETURN> for help text, or \
mayavi_utils( for pop-up help in Canopy")
 
@ml.animate(delay=30)
def spin():
    '''
    Spin the view interactively
    ''' 
    f = ml.gcf()
    while 1:
        f.scene.camera.azimuth(1)
        f.scene.render()
        yield

@ml.animate(delay=10,ui=None)
def zoom(e=1.0,n=200):
    '''
    Zoom the view a factor exp(e), in n steps
    ''' 
    f=np.exp(e)**(1.0/n)
    g=ml.gcf()
    for i in range(n):
        g.scene.camera.zoom(f)
        g.scene.render()
        yield

@ml.animate(delay=10)
def spin_zoom(factor=2.,n=360,phi=360,name='im',save=0):
    '''
    Spin + zoom the view interactively
    ''' 
    f = ml.gcf()
    eps=factor/n
    dphi=phi/n
    for i in range(n):
        f.scene.camera.azimuth(dphi)
        f.scene.camera.zoom(1.0+eps)
        f.scene.render()
        if save:
            file=name+'_{:04d}.png'.format(i)
            ml.savefig(file)
        yield
    for i in range(n):
        f.scene.camera.azimuth(dphi)
        f.scene.camera.zoom(1.0-eps)
        f.scene.render()
        if save:
            file=name+'_{:04d}.png'.format(i+n)
            ml.savefig(file)
        yield

@ml.animate(delay=10,ui=None)
def dolly(e=1.0,n=200):
    '''
    Dolly the view a factor exp(e), in n steps
    ''' 
    f=np.exp(e)**(1.0/n)
    g=ml.gcf()
    for i in range(n):
        g.scene.camera.dolly(f)
        g.scene.render()
        yield

@ml.animate(delay=10)
def spin_movie(dir='anim',n=360):
    '''
    Spin the view 360 degrees, in n steps, saving images in dir
    ''' 
    os.system('mkdir -p '+dir)
    g = ml.gcf()
    dphi=360.0/n
    for i in range(n):
        g.scene.camera.azimuth(dphi)
        g.scene.render()
        ml.show()
        file=dir+'/im_{:03d}.png'.format(i)
        print(file)
        yield
        #ml.savefig(file)

def cube(x=0,y=0,z=0,size=1,phi=30):
    '''
    Construct a 3D cube
    '''
    c=np.cos(phi+np.pi/4)
    s=np.sin(phi+np.pi/4)
    x1=0.5*size*c*np.sqrt(2.)
    y1=0.5*size*s*np.sqrt(2.)
    z1=0.5*size
    x2=x+np.array([+x1,-y1,-x1,+y1,+x1])
    y2=y+np.array([+y1,+x1,-y1,-x1,+y1])
    z2=z+np.array([+z1,+z1,+z1,+z1,+z1])
    lw=0.1
    ml.plot3d(x2,y2,+z2,tube_radius=None,line_width=lw)
    ml.plot3d(x2,y2,-z2,tube_radius=None,line_width=lw)
    for i in range(4):
        ml.plot3d([x2[i],x2[i]],[y2[i],y2[i]],[-z1,z1],tube_radius=None,line_width=lw)

def torus(size=1.,levels=1,nphi=31,n=31,azimuth=90.):
    '''
    Construct a torus consisting of cubes with given size
    '''
    dphi=2.*np.pi/nphi
    if levels>1:
        size=2.2*0.9**(1./3.)*3**levels
    radius=size/dphi
    for i in range(n):
        phi=i*dphi+np.pi/2.
        x=radius*np.cos(phi)
        y=radius*(np.sin(phi)-1.0)
        z=0.0
        cube(x,y,z,size,phi)        
    ml.draw()
    ml.view(focalpoint=[0.,0.,0.],azimuth=azimuth,elevation=75., 
      distance=4*radius)
