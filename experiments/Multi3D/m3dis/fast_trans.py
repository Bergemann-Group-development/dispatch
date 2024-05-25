from .trans import trans_class
import numpy as np
import pandas as pd

class fast_trans(trans_class):
  def __init__(self, kr, atom):
    self.kr = kr
    self.atom = atom
    self.el = atom
    self.bbt = atom.bbt
    self.title = 'Line = {:.2f} [Å]'.format(self.lam0)
  
  def __getattr__(self, attr):
    if attr in self.bbt:
      return self.bbt.at[self.kr, attr]
    else:
      raise AttributeError('fast_trans object does not have attribute ', attr)
    
  def __repr__(self):
    return self.bbt.loc[self.kr].__repr__()