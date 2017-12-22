'''
   Staggered derivatives and interpolations, corresponding to the ones used in the stagger2
   solver in the DISPATCH code
'''
from numpy import roll, log, exp

def xdn(f):
    return (f+roll(f,1,0))*0.5
def ydn(f):
    return (f+roll(f,1,1))*0.5
def zdn(f):
    return (f+roll(f,1,2))*0.5
def xup(f):
    return (f+roll(f,-1,0))*0.5
def yup(f):
    return (f+roll(f,-1,1))*0.5
def zup(f):
    return (f+roll(f,-1,2))*0.5

def ddxdn(f,s):
    c=1./s.ds[0]
    return (f-roll(f,1,0))*c
def ddydn(f,s):
    c=1./s.ds[1]
    return (f-roll(f,1,1))*c
def ddzdn(f,s):
    c=1./s.ds[2]
    return (f-roll(f,1,2))*c

def ddxup(f,s):
    c=1./s.ds[0]
    return (roll(f,-1,0)-f)*c
def ddyup(f,s):
    c=1./s.ds[1]
    return (roll(f,-1,1)-f)*c
def ddzup(f,s):
    c=1./s.ds[2]
    return (roll(f,-1,2)-f)*c

def ddx(f,s):
    c=2./s.ds[0]
    return (roll(f,-1,0)-roll(f,1,0))*c
def ddy(f,s):
    c=2./s.ds[1]
    return (roll(f,-1,1)-roll(f,1,1))*c
def ddz(f,s):
    c=2./s.ds[2]
    return (roll(f,-1,2)-roll(f,1,2))*c

# These functions may be used to calculate the velocity from the momenta and mass density
def fux(f):
    return xup(f.px/exp(xdn(log(f.d))))
def fuy(f):
    return yup(f.py/exp(ydn(log(f.d))))
def fuz(f):
    return zup(f.pz/exp(zdn(log(f.d))))
