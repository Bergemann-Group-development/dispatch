'''
   compact data reader
'''
import numpy as np
import read_file as rf
reload(rf)
from read_file import read_file

def read_data (ip=1,it=0,level=0,ix=0,iy=0,iz=0,ng=0,st=0,nodata=0,run='',dir='../data'):
    if level>0:
        ip=14+ix+3*iy+9*iz+27*(level-1)
    name=dir+'/'+run+'/{:04d}_{:04d}'.format(ip,it)
    f=read_file(name,nodata=nodata)
    if not f:
        name=dir+'/'+run+'/{:04d}_{:04d}'.format(ip,it)
        f=read_file(name,nodata=nodata)
    n=f.n
    f.x=f.pos[0]+f.ds[0]*(np.arange(f.n[0])-0.5*(f.n[0]-1.0))
    f.y=f.pos[1]+f.ds[1]*(np.arange(f.n[1])-0.5*(f.n[1]-1.0))
    f.z=f.pos[2]+f.ds[2]*(np.arange(f.n[2])-0.5*(f.n[2]-1.0))
    return f
