import numpy as np
import matplotlib.pyplot as pl
import fortranfile as ftf
import dispatch_utils as du
import scaling
from time import time

class particle():
    def __init__(self,a):
        self.p=a[0:3]
        self.t=a[6]
        self.s=a[7]
        self.w=a[8]

def mass(it,run,data):
    ''' compute mass and particle dict '''
    print('snapshot',it)
    files=du.patch_files(iout=it,run=run,data=data)
    dict={}
    t0=time()
    m=0.0
    npe=0
    for file in files:
        f=ftf.FortranFile(file+'.peb','rb')
        fmt=f.readInts()
        dim=f.readInts()
        a=f.readReals().reshape(dim[1],dim[0]).transpose()
        idv=f.readInts()
        f.close()
        idx=idv[0:dim[0]]
        npe+=dim[0]
        for i in range(dim[0]):
            dict[idx[i]]=a[i,8]
            m+=a[i,8]
    print('snapshot',it,m,time()-t0,npe)
    return m,dict

def delta(it,run,data,dict):
    ''' accretion relative to previous particle dict '''
    files=du.patch_files(iout=it,run=run,data=data)
    t0=time()
    t2=t0
    m=0.0
    ma=0.0
    mo=0.0
    npe=0
    no=0
    na=0
    mindist=1e9
    accr={}
    for file in files:
        f=ftf.FortranFile(file+'.peb','rb')
        fmt=f.readInts()
        dim=f.readInts()
        a=f.readReals().reshape(dim[1],dim[0]).transpose()
        idv=f.readInts()
        f.close()
        idx=idv[0:dim[0]]
        npe+=dim[0]
        
        # loop over all indices
        for i in range(dim[0]):
            try:
                m+=dict[idx[i]]
            
            # if doesn't exist in next snapshot it has been accreted or lost
            except:
                dist=np.sqrt(np.sum(a[i,0:3]**2))
                # accretion if from small distances
                if dist < 200:
                    ma+=a[i,8]
                    na+=1
                    dict[idx[i]]=a[i,8]
                    accr[idx[i]]=particle(a[i,:])
                # lost at large distances
                else:
                    no+=1
                    mo+=a[i,8]
                    if dist<mindist:
                        mindist=dist
                        #print 'min',mindist
    print 'snapshot',it,m,ma,mo,time()-t2,npe,no,na,mindist
    return ma,dict,accr

def scan(it0,it1,run,data):
    ''' scan accretion over a range of snapshots '''
    m,pebs=mass(it1,run,data)
    dict={}
    for it in range(it1-1,it0-1,-1):
        dm,pebs,accr=delta(it,run,data,pebs)
        dict[it]=accr
        s=cgs.r_earth*np.array([p.s for p in accr.values()])
        pl.hist(np.log10(s))
        pl.title(it)
    return dict
