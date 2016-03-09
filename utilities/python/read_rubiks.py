'''
   Compact data reader for Rubik's cube arrangement, which returns both
   the data and appropriately shifted and scaled meshes.  The position
   may be either given as a patch number, or via ix,iy,iz,level combo
'''
import numpy as np
from read_file import read_file

def read_rubiks (ip=1,it=0,level=0,ix=0,iy=0,iz=0,ng=0,run='pebble_atmosphere',dir='data'):
    if level>0:
        ip=14+ix+3*iy+9*iz+27*(level-1)
    name=dir+'/'+run+'_{:03d}_{:04d}'.format(ip,it)
    f=read_file(name)
    if not f:
        name=dir+'/'+run+'_{:04d}_{:04d}'.format(ip,it)
        f=read_file(name)
    n=f.d.shape
    def mesh(id,ng=0):
        print "id =",id, "ng =",ng
        level=(id-1)/27
        i=id-level*27-1
        class m:
            level=(id-1)/27+1
            ix=np.mod(i,3)-1
            iy=np.mod(i/3,3)-1
            iz=i/9-1
            s=3.0**(1-level)
            x=(np.linspace(-s/2,s/2,n[0])+ix*s)*(n[0]-1.0+2.*ng)/n[0]
            y=(np.linspace(-s/2,s/2,n[1])+iy*s)*(n[1]-1.0+2.*ng)/n[1]
            z=(np.linspace(-s/2,s/2,n[2])+iz*s)*(n[2]-1.0+2.*ng)/n[2]
        return m
    m=mesh(ip,ng=ng)
    f.x=m.x
    f.y=m.y
    f.z=m.z
    return f
