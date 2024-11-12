from .trans import trans_class
import numpy as np
import pandas as pd

class fast_trans(trans_class):
  #=====================================================================
  def __init__(self, kr, atom):
    #-----------------------------------------------------------------
    self.kr = kr
    self.atom = atom
    self.bbt = atom.bbt
    self.title = 'Line = {:.2f} [Å]'.format(self.lam0)
    self.prepared = False
  
  def __getattr__(self, attr):
    if attr in self.bbt:
      return self.bbt.at[self.kr, attr]
    else:
      raise AttributeError('fast_trans object does not have attribute ', attr)
    
  def __repr__(self):
    return self.bbt.loc[self.kr].__repr__()

  def prepare(self):
    if hasattr(self, 'prepared'):
      if self.prepared: return
      
    self.run = self.atom.run
  
class llist_trans(fast_trans):
  #=====================================================================
  def __init__(self, kr, run):
    #-----------------------------------------------------------------
    self.kr = kr
    self.run = run
    self.bbt = run.llist
    self.title = 'Line = {:.2f} [Å], {:}'.format(self.lam0, self.ion)
    self.prepared = False
    self.lam_blue = self.lam_vac - 1.0
    self.lam_red  = self.lam_vac + 1.0

  def prepare(self):
    return trans_class.prepare(self)