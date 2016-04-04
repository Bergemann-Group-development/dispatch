'''
   compact data reader
'''
import numpy as np
import read_file as rf
reload(rf)
from read_file import read_file

def read_data (ip=1,it=0,level=0,ix=0,iy=0,iz=0,ng=0,st=0,run='pebble_atmosphere',dir='data'):
    if level>0:
        ip=14+ix+3*iy+9*iz+27*(level-1)
    name=dir+'/'+run+'/{:03d}_{:04d}'.format(ip,it)
    f=read_file(name)
    if not f:
        name=dir+'/'+run+'/{:04d}_{:04d}'.format(ip,it)
        f=read_file(name)
    n=f.d.shape
    def mesh(id,ng=0):
        #print "id =",id, "ng =",ng
        level=(id-1)/27
        i=id-level*27-1
        class m:
            level=(id-1)/27+1
            ix=np.mod(i,3)-1
            iy=np.mod(i/3,3)-1
            iz=i/9-1
            s=3.0**(1-level)
            x=(np.linspace(-s/2,s/2,n[0]))*(n[0]-1.0)/(n[0]-2*ng) + ix*s
            y=(np.linspace(-s/2,s/2,n[1]))*(n[1]-1.0)/(n[1]-2*ng) + iy*s
            z=(np.linspace(-s/2,s/2,n[2]))*(n[2]-1.0)/(n[2]-2*ng) + iz*s
            if st==1: x=x-0.5*s/n[0]
            if st==2: y=y-0.5*s/n[1]
            if st==3: z=z-0.5*s/n[2]
        return m
    m=mesh(ip,ng=ng)
    f.x=m.x
    f.y=m.y
    f.z=m.z
    return f