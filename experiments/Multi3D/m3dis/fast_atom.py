import io, types
import pandas as pd
from .auxfuncs import *

def read(text, dtype):
  return np.array(next(text).strip().split()).astype(dtype)

def read_table(file, **kwargs):
  return pd.read_csv(file, sep='\s+', header=None, comment='*', quotechar="'", **kwargs)


class fast_atom(object):
  def __init__(self, atom_file, qnorm):
    f = open(atom_file, 'r')
    condition = lambda x: not x.strip().startswith('*')
    text = filter(condition, f)
    self.id = next(text).strip()
    self.molec = self.id in ['H2', 'C2', 'N2', 'O2', 'CH', 'OH', 'NH', 'CN', 'CO', 'NO', 'H2O']
    if not self.molec: self.id = self.id.capitalize()

    self.abnd, self.ai = read(text, float)
    self.awgt = self.ai * UU
    self.nk, self.nline, self.ncont, self.nrfix = read(text, int)
    self.nrad = self.nline + self.ncont
    f.close()
    levels = read_table(atom_file, skiprows=3, nrows=self.nk, names=['ecm', 'g', 'label', 'ion'])
    bb     = read_table(atom_file, skiprows=3+self.nk, nrows=self.nline, names=['j', 'i', 'f', 'nq', 'qmax', 'q0', 'iw', 'ga', 'gw', 'gq', 'profile'])
    levels['ev'] = levels.ecm * CC*HH/EE
    bb.i = bb.i-1
    bb.j = bb.j-1
    Ei = levels.iloc[bb.i].ecm.values
    Ej = levels.iloc[bb.j].ecm.values
    bb['kr'] = range(self.nline)
    bb['nu0'] = CC * (Ej - Ei)
    bb['lam_vac'] = 1e8 * CC / bb.nu0
    bb['lam0'] = vacuum2obs(bb.lam_vac)
    bb['loggf'] = np.log10(bb.f*levels.iloc[bb.i].g.values)
    self.qnorm = qnorm * 1e5 / CC
    qmax = get_qmax(bb.qmax, bb.q0) * self.qnorm
    bb['lam_vac_blue'] = bb.lam_vac / (1 + qmax)
    bb['lam_vac_red']  = bb.lam_vac / (1 - qmax)
    bb['lb'] = vacuum2obs(bb.lam_vac_blue)
    bb['lr'] = vacuum2obs(bb.lam_vac_red)

    self.level = levels
    self.bb = bb