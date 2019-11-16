"""
    Utility to read a single data/run/SSSSS/aux_PPPPP.dat file.  Return
    
    aux.filename        ! filename
    aux.version!        ! aux file version
    aux.vars            ! dict with variable_names,variable_values
"""
from scipy.io import FortranFile
import numpy as np
import os

def to_str(ii):
    s=''
    for i in ii:
        if i != 32:
            s+=str(chr(i))
    return s

class aux():
  """ aux format file object """
  def __init__(self,id=1,io=0,run='',data='../data',file=None,verbose=0):
    if os.path.isfile(file):
        self.filename=file
    else:
        self.filename=data+'/'+run+'/{:05d}/{:05d}.aux'.format(io,id)
    self.verbose=verbose
    self._read()

  def _read(self):
    """ read an aux format file """
    with FortranFile(self.filename,'r') as ff:
      self.version,self.id=ff.read_ints()                # version, id
      if self.verbose>0:
        print('version {}'.format(self.version))
        print('patch: {}'.format(self.id))
      vars={}
      while True:
        try:
          ii=ff.read_record('b')                         # name
          name=to_str(ii)
          if self.verbose>1:
            print(' name:',name)
          rnk=ff.read_ints()[0]                          # shape
          if self.verbose>2:
            print('rank:',rnk)
          shp=ff.read_ints()
          if self.verbose>2:
            print('shape:',shp)
          rank=len(shp)
          ii=ff.read_record('b')                         # type
          tpe=to_str(ii)
          if self.verbose>2:
            print(' type:',tpe)
          if tpe[0]=='r':
            v=ff.read_reals('<f4')                       # value
          else:
            v=ff.read_ints('<i4')
          v=v.reshape(shp[-1::-1]).transpose()
          if self.verbose>3:
            print('  min: {}\n  max: {}'.format(v.min(),v.max()))
          vars[name]={'type':tpe,'name':name,'rank':rank,'shape':shp,'v':v}
        except:
          break
    self.vars=vars

  def var(self,k):
    return self.vars[k]['v']

if __name__=='__main__':
    with FortranFile('test.aux','w') as ff:
        ff.write_record((1,3))

        ff.write_record('test1')
        shp=(2,3)
        ff.write_record(shp)
        ff.write_record('r3')
        v=np.ones(shp,dtype=np.float32)
        ff.write_record(v)

        ff.write_record('test2')
        shp=(2,4,5)
        ff.write_record(shp)
        ff.write_record('i3')
        v=np.ones(shp,dtype=np.int32)
        ff.write_record(v)
    test=aux('test.aux',verbose=4)
