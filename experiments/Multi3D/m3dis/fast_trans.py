from .trans import trans_class
import numpy as np
import pandas as pd

class fast_trans(trans_class, pd.core.series.Series):
  def __init__(self, Series):
    self.__dict__ = Series.__dict__

class lines_class(object):
  def __init__(self, atom, df, run):
    df['i0'] = np.searchsorted(run.lam, df.lam0)
    df['ib'] = np.searchsorted(run.lam, df.lb)
    df['ir'] = np.searchsorted(run.lam, df.lr)

    if hasattr(run, 'lvlerr'):
        df['ilvlerr'] = run.lvlerr[df.i]
        df['jlvlerr'] = run.lvlerr[df.j]

    self.df = df
    self.run = run
    self.atom = atom

    self.Series = [None,] * len(df)
    #if prep: self.Series = [fast_trans(i) for i in df.iloc]

  def __len__(self):
    return len(self.df)
  
  def __repr__(self):
    return self.df.__repr__()
  
  def _repr_html_(self):
    return self.df._repr_html_()
  
  def __getitem__(self, index):
    #if hasattr(self, 'Series'):
    #  ser = self.Series[index]
    #else:
    #  ser = fast_trans(self.df.iloc[index])
    if self.Series[index] is None:
      self.Series[index] = fast_trans(self.df.iloc[index])

    ser = self.Series[index]
    ser.title = 'Line = {:.2f} [Å]'.format(ser.lam0)
    ser.run = self.run
    ser.el = self.atom
    return ser