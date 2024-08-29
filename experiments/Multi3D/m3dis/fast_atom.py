import io, types
import pandas as pd
from .fast_trans import fast_trans
from .auxfuncs import *
from .atom import atom_class

formatters = {'ga':'{:.3e}'.format,
              'gv':'{:.3e}'.format,
              'gq':'{:.2e}'.format,
              'ev':'{:.6f}'.format,
              'evi':'{:.6f}'.format,
              'lam_air':'{:.3f}'.format,
              'alpha':'{:.3e}'.format,
              'label':"'{:}'".format,
              'loggf':'{:.3f}'.format}

def read(text, dtype):
  return np.array(next(text).strip().split()).astype(dtype)

def read_table(file, **kwargs):
  return pd.read_csv(file, sep=r'\s+', comment='*', quotechar="'", **kwargs)

class fast_atom(object):
  def __init__(self, atom_file):
    with open(atom_file, 'r') as f:
      condition = lambda x: not x.strip().startswith('*')
      text = filter(condition, f)
      self.id = next(text).strip()
      self.molec = self.id in ['H2', 'C2', 'N2', 'O2', 'CH', 'OH', 'NH', 'CN', 'CO', 'NO', 'H2O']
      if not self.molec: self.id = self.id.capitalize()

      self.abnd, self.ai = read(text, float)
      self.awgt = self.ai * UU
      self.nk, self.nline, self.ncont, self.nrfix = read(text, int)
      self.nrad = self.nline + self.ncont
      
    levels = read_table(atom_file, header=2, nrows=self.nk, names=['ecm', 'g', 'label', 'ion'])
    bb     = read_table(atom_file, header=2+self.nk, nrows=self.nline, names=['j', 'i', 'f', 'nq', 'qmax', 'q0', 'iw', 'ga', 'gv', 'gq', 'profile'])
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
    pass

  def post_init(self):
    self.molec = self.id in ['H2', 'C2', 'N2', 'O2', 'CH', 'OH', 'NH', 'CN', 'CO', 'NO', 'H2O']
    if not self.molec: self.id = self.id.capitalize()
    self.awgt = self.ai * UU

    self.lvl_tab = self.lvl.columns
    if hasattr(self, 'bft'): self.bft_tab = self.bft.columns
    if hasattr(self, 'col'): self.col_tab = self.col.columns
    if hasattr(self, 'bbt'):
      self.bbt_tab = self.bbt.columns
      init = np.frompyfunc(fast_trans, 2, 1)
      self.trans = init(np.arange(self.nline), self)
  
  @classmethod
  def read(cls, atom_file):
    self = cls()
    
    with open(atom_file, 'r') as f:
      condition = lambda x: not x.strip().startswith('*')
      text = filter(condition, f)
      self.id = next(text).strip()
      self.abnd, self.ai = read(text, float)
      self.nk, self.nline, nnu_bf, self.ncol = read(text, int)

    pdkwargs = {'filepath_or_buffer':atom_file, 'sep':r'\s+', 'quotechar':"'", 'comment':'*'}
    pdkwargs['header'] = 3
    self.lvl = pd.read_csv(nrows=self.nk, **pdkwargs)
    pdkwargs['header'] = 3+self.nk+1
    self.bbt = pd.read_csv(nrows=self.nline, **pdkwargs)
    pdkwargs['header'] = 3+self.nk+1+self.nline+1
    if nnu_bf > 0:
      self.bft = pd.read_csv(nrows=nnu_bf, **pdkwargs)
    if self.ncol > 0:
      pdkwargs['header'] = 3+self.nk+1+self.nline+1+nnu_bf+1
      self.col = pd.read_csv(nrows=self.ncol, **pdkwargs)
    self.bbt = self.bbt.rename(columns={'lam_air':'lam0'})
    self.post_init()
    return self
  
  @classmethod
  def from_old(cls, atom_file):
    self = cls()
    atom = atom_class(atom_file, qnorm=10., read_bf=True, read_coll=True, sanity_check=True)
    fast = fast_atom(atom_file)
    self.id = atom.id
    self.abnd = atom.abnd
    self.ai = atom.ai
    self.nk = atom.nk
    self.nline = atom.nline
    self.lvl = fast.level[['ev', 'g', 'ion', 'label']]

    bb = fast.bbt.copy()
    bb['evi'] = self.lvl.ev[bb.i].values
    bb = bb[['i', 'j', 'lam0', 'evi', 'loggf', 'ga', 'gv', 'gq', 'nq', 'qmax', 'q0', 'profile']]
    bb.i = bb.i
    bb.j = bb.j
    self.bbt = bb

    bf = pd.DataFrame()
    cont = atom.trans[atom.nline:]
    if len(cont) > 0:
      bf['i'] = np.concatenate([[c.i for i in range(c.lam.size)] for c in cont])
      bf['lam_air'] = np.concatenate([c.lam for c in cont])
      bf['alpha'] = np.concatenate([c.alfac for c in cont])
      self.bft = bf

    if len(atom.colls) > 0:
      ctemp = atom.colls[0].coltemp
      crate = np.array([c.colrate for c in atom.colls])
      cl = pd.DataFrame(data=crate, columns=ctemp)
      cl['key'] = [c.key for c in atom.colls]
      cl['i'] = [c.lower.i for c in atom.colls]
      cl['j'] = [c.upper.i for c in atom.colls]
      index = list(cl.columns[-3:]) + list(cl.columns[:-3])
      self.col = cl[index]

    self.post_init()
    return self
  
  @classmethod
  def from_pandas(cls, id, abnd, ai, lvl, bbt=None, bft=None, col=None):
    self = cls()
    self.id = id
    self.abnd = abnd
    self.ai = ai
    self.lvl = lvl
    self.nk = len(lvl)
    if bbt is not None:
      self.bbt = bbt
      self.nline = len(bbt)
    if bft is not None:
      self.bft = bft
    if col is not None:
      self.col = col
      self.ncol = len(col)
    self.post_init()
    return self

  def write(self, out_file):
    lvl = self.lvl[self.lvl_tab].sort_values('ev').reset_index()
    index_map = pd.Series(data=lvl.index, index=lvl['index'])
    lvl = lvl.drop(columns=['index'])

    if hasattr(self, 'bbt'):
      bbt = self.bbt[self.bbt_tab].sort_values('lam0', ignore_index=True)
      bbt['i'] = index_map[bbt.i].values
      bbt['j'] = index_map[bbt.j].values
      bbt = bbt.rename(columns={'lam0':'lam_air'})
    else:
      bbt = []

    if hasattr(self, 'bft'):
      bft = self.bft[self.bft_tab]
      bft['i'] = index_map[bft.i].values
      bft = bft.sort_values(['i', 'lam_air'], ascending=[True, False], ignore_index=True)
    else:
      bft = []

    if hasattr(self, 'col'):
      for t in self.col.columns[3:].values:
        formatters[t]='{:.4e}'.format

      col = self.col[self.col_tab]
      col['i'] = index_map[col.i].values
      col['j'] = index_map[col.j].values
      col = col.sort_values(['key', 'i', 'j'], ascending=[True, True, True], ignore_index=True)
    else:
      col = []

    f = open(out_file, 'w+')
    f.write(f'{self.id:}\n')
    f.write(f'{self.abnd:} {self.ai:}\n')
    f.write(f'{len(lvl):} {len(bbt):} {len(bft):} {len(col):}\n')
    f.write(lvl.to_string(formatters=formatters) + '\n')
    if len(bbt)>0: f.write(bbt.to_string(formatters=formatters) + '\n')
    if len(bft)>0: f.write(bft.to_string(formatters=formatters, index=False) + '\n')
    if len(col)>0: f.write(col.to_string(formatters=formatters) + '\n')
    f.close()

  def load_sp(self, lam, qnorm):
    self.qnorm = qnorm * 1e5 / CC
    qmax = get_qmax(self.bbt.qmax, self.bbt.q0) * self.qnorm
    self.bbt['lb'] = self.bbt.lam0 / (1 + qmax)
    self.bbt['lr'] = self.bbt.lam0 / (1 - qmax)
    self.bbt['i0'] = np.searchsorted(lam, self.bbt.lam0)
    self.bbt['ib'] = np.searchsorted(lam, self.bbt.lb)
    self.bbt['ir'] = np.searchsorted(lam, self.bbt.lr)

  def to_linelist(self, file):
    f = open(file, 'w+')
    ions = np.unique(self.lvl.ion)
    for ion in ions:
      ll = pd.DataFrame()
      mask = np.array(self.lvl.ion[self.bbt.i] == ion)
      ev0 = self.lvl.ev[self.lvl.ion == ion].values[0]
      if any(mask):
        f.write(f"'  {np.round(self.ai).astype(int):<16.3f}' {ion:d} {np.sum(mask):d}\n")
        f.write(f"'{self.id:3} {'I'*ion:<3}'\n")
        bbt = self.bbt[mask]
        ll['lam0'] = bbt['lam0']
        ll['evi'] = bbt['evi'] - ev0
        ll['loggf'] = bbt['loggf']
        ll['g_j'] = self.lvl.g[self.bbt.j[mask]].values
        ll['gv'] = bbt['gv']
        ll['ga'] = bbt['ga']
        ll['gq'] = bbt['gq']
        ll['lower_orb'] = 'x'
        ll['upper_orb'] = 'x'
        ll['ent8'] = 0
        ll['ent9'] = 0
        ll['label'] = 'x'
        f.write(ll.to_string(index=False, header=False))

    f.close()

  def remove_levels(self, ilvl):
    dropped = self.lvl.loc[ilvl]
    self.lvl = self.lvl.drop(ilvl).reset_index()
    index_map = pd.Series(data=self.lvl.index, index=self.lvl['index'])
    self.lvl = self.lvl.drop(columns=['index'])
    self.nlvl = self.lvl.shape[0]
    if hasattr(self, 'bbt'): self.bbt['i'] = index_map[self.bbt.i].values
    if hasattr(self, 'bbt'): self.bbt['j'] = index_map[self.bbt.j].values
    if hasattr(self, 'bft'): self.bft['i'] = index_map[self.bft.i].values
    if hasattr(self, 'col'): self.col['i'] = index_map[self.col.i].values
    if hasattr(self, 'col'): self.col['j'] = index_map[self.col.j].values
    return dropped