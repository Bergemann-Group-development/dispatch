from .trans import Transition
import numpy as np
import pandas as pd
import functools

class TransFactory(object):
  def __init__(self, atom, in_lines=None):
    self.atom = atom
    self.in_lines = in_lines or slice(None, None)
    self.krs = np.arange(len(atom.bbt))[self.in_lines]

  @functools.cache
  def __getitem__(self, key):
    kr = self.krs[key]
    return FastTransition(kr, self.atom)
  
  def get_lines(self, in_lines):
    return TransFactory(self.atom, in_lines=in_lines)


class FastTransition(Transition):
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
      raise AttributeError('FastTransition object does not have attribute ', attr)
    
  def __repr__(self):
    return self.bbt.loc[self.kr].__repr__()

  def prepare(self):
    if hasattr(self, 'prepared'):
      if self.prepared: return
      
    self.run = self.atom.run