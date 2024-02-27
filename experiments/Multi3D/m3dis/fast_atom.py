import io, types
import pandas as pd
from .fast_trans import fast_trans
from .auxfuncs import *

formatters = {'ga':'{:.3e}'.format,
              'ev':'{:.6f}'.format,
              'gq':'{:.2e}'.format,
              'lam_air':'{:.3f}'.format,
              'alpha':'{:.3e}'.format,
              'label':"'{:}'".format,
              'loggf':'{:.3f}'.format}

def read(text, dtype):
  return np.array(next(text).strip().split()).astype(dtype)

def read_table(file, **kwargs):
  return pd.read_csv(file, sep='\s+', comment='*', quotechar="'", **kwargs)

class fast_atom(object):
  def __init__(self, atom_file):
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
    levels = read_table(atom_file, header=2, nrows=self.nk, names=['ecm', 'g', 'label', 'ion'])
    bb     = read_table(atom_file, header=2+self.nk, nrows=self.nline, names=['j', 'i', 'f', 'nq', 'qmax', 'q0', 'iw', 'ga', 'gw', 'gq', 'profile'])
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

    self.level = levels
    self.bbt = bb

    init = np.frompyfunc(fast_trans, 2, 1)
    self.trans = init(np.arange(self.nline), self)

  def load_sp(self, lam, qnorm):
    self.qnorm = qnorm * 1e5 / CC
    qmax = get_qmax(self.bbt.qmax, self.bbt.q0) * self.qnorm
    self.bbt['lb'] = self.bbt.lam0 / (1 + qmax)
    self.bbt['lr'] = self.bbt.lam0 / (1 - qmax)
    self.bbt['i0'] = np.searchsorted(lam, self.bbt.lam0)
    self.bbt['ib'] = np.searchsorted(lam, self.bbt.lb)
    self.bbt['ir'] = np.searchsorted(lam, self.bbt.lr)

class pandas_atom(object):
  def __init__(self):
    self.molec = self.id in ['H2', 'C2', 'N2', 'O2', 'CH', 'OH', 'NH', 'CN', 'CO', 'NO', 'H2O']
    if not self.molec: self.id = self.id.capitalize()
    self.awgt = self.ai * UU

    self.lvl_tab = self.lvl.columns
    self.bbt_tab = self.bbt.columns
    self.bft_tab = self.bft.columns
    self.col_tab = self.col.columns

    init = np.frompyfunc(fast_trans, 2, 1)
    self.trans = init(np.arange(self.nline), self)
  
  @classmethod
  def read(cls, atom_file):
    f = open(atom_file, 'r')
    cls.id = next(f).strip()
    cls.abnd, cls.ai = read(f, float)
    cls.nk, cls.nline, nnu_bf, cls.ncol = read(f, int)
    f.close()

    pdkwargs = {'filepath_or_buffer':atom_file, 'sep':'\s+', 'quotechar':"'", 'comment':'*'}
    pdkwargs['header'] = 3
    cls.lvl = pd.read_csv(nrows=cls.nk, **pdkwargs)
    pdkwargs['header'] = 3+cls.nk+1
    cls.bbt = pd.read_csv(nrows=cls.nline, **pdkwargs)
    pdkwargs['header'] = 3+cls.nk+1+cls.nline+1
    cls.bft = pd.read_csv(nrows=nnu_bf, **pdkwargs)
    pdkwargs['header'] = 3+cls.nk+1+cls.nline+1+nnu_bf+1
    cls.col = pd.read_csv(nrows=cls.ncol, **pdkwargs)
    cls.bbt = cls.bbt.rename(columns={'lam_air':'lam0'})
    return cls()
  
  @classmethod
  def from_old(cls, atom_file):
    atom = atom_class(atom_file, qnorm=10., read_bf=True, read_coll=True, sanity_check=True)
    fast = fast_atom(atom_file, qnorm=10.)
    cls.id = atom.id
    cls.abnd = atom.abnd
    cls.ai = atom.ai
    cls.nk = atom.nk
    cls.lvl = fast.level[['ev', 'g', 'ion', 'label']]

    bb = fast.bbt.copy()
    bb = bb[['lam0', 'i', 'j', 'loggf', 'nq', 'qmax', 'q0', 'ga', 'gw', 'gq', 'profile']]
    bb.i = bb.i+1
    bb.j = bb.j+1
    cls.bbt = bb

    bf = pd.DataFrame()
    cont = atom.trans[atom.nline:]
    bf['i'] = np.concatenate([[c.i+1 for i in range(c.lam.size)] for c in cont])
    bf['lam_air'] = np.concatenate([c.lam for c in cont])
    bf['alpha'] = np.concatenate([c.alfac for c in cont])
    cls.bft = bf

    ctemp = atom.colls[0].coltemp
    crate = np.array([c.colrate for c in atom.colls])
    cl = pd.DataFrame(data=crate, columns=ctemp)
    cl['key'] = [c.key for c in atom.colls]
    cl['i'] = [c.lower.i+1 for c in atom.colls]
    cl['j'] = [c.upper.i+1 for c in atom.colls]
    index = list(cl.columns[-3:]) + list(cl.columns[:-3])
    cls.col = cl[index]
    return cls()

  def write(self, out_file):
    for t in self.col.columns[3:].values:
      formatters[t]='{:.2e}'.format
    lvl = self.lvl[self.lvl_tab].sort_values('ev', ignore_index=True)
    bbt = self.bbt[self.bbt_tab].sort_values('lam0', ignore_index=True)
    bft = self.bft[self.bft_tab].sort_values(['i', 'lam_air'], ascending=[True, False], ignore_index=True)
    col = self.col[self.col_tab].sort_values(['key', 'i', 'j'], ascending=[True, True, True], ignore_index=True)

    bbt = bbt.rename(columns={'lam0':'lam_air'})

    f = open(out_file, 'w+')
    f.write(f'{self.id:}\n')
    f.write(f'{self.abnd:} {self.ai:}\n')
    f.write(f'{len(lvl):} {len(bbt):} {len(bft):} {len(col):}\n')
    f.write(lvl.to_string(formatters=formatters) + '\n')
    f.write(bbt.to_string(formatters=formatters) + '\n')
    f.write(bft.to_string(formatters=formatters, index=False) + '\n')
    f.write(col.to_string(formatters=formatters) + '\n')
    f.close()

  def load_sp(self, lam, qnorm):
    self.qnorm = qnorm * 1e5 / CC
    qmax = get_qmax(self.bbt.qmax, self.bbt.q0) * self.qnorm
    self.bbt['lb'] = self.bbt.lam0 / (1 + qmax)
    self.bbt['lr'] = self.bbt.lam0 / (1 - qmax)
    self.bbt['i0'] = np.searchsorted(lam, self.bbt.lam0)
    self.bbt['ib'] = np.searchsorted(lam, self.bbt.lb)
    self.bbt['ir'] = np.searchsorted(lam, self.bbt.lr)