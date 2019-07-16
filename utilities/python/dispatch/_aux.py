"""
    Utility to read a single data/run/SSSSS/aux_PPPPP.dat file.  Return
    
    aux.filename        ! filename
    aux.version!        ! aux file version
    aux.vars            ! dict with variable_names,variable_values
"""
from scipy.io import FortranFile
import numpy as np

class aux():
  """ aux format file object """
  def __init__(self,filename,verbose=0):
    self.filename=filename
    self.verbose=verbose
    self._vars()

  def _vars(self):
    """ read an aux format file """
    with FortranFile(self.filename,'r') as ff:
      self.version,self.id=ff.read_ints()
      if self.verbose>0:
        print('version={}'.format(self.version))
      vars={}
      while True:
        try:
          name=ff.read_record('U5')[0]
          if self.verbose>1:
            print(' name:',name)
          shp=ff.read_ints()
          if self.verbose>2:
            print('shape:',shp)
          rank=len(shp)
          tpe=ff.read_record('U2')[0]
          if self.verbose>2:
            print(' type:',tpe)
          if tpe[0]=='r':
            v=ff.read_reals('<f4')
          else:
            v=ff.read_ints('<i4')
          v=v.reshape(shp)
          if self.verbose>3:
            print('  min: {}\n  max: {}'.format(v.min(),v.max()))
          vars[name]={'type':tpe,'name':name,'rank':rank,'shape':shp,'v':v}
        except:
          break
    self.vars=vars

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