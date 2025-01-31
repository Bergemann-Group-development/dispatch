from m3dis.auxfuncs import *
import numpy as np
import pandas as pd

element_ids = [
        'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F ',
        'Ne','Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar',
        'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co',
        'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
        'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh',
        'Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',
        'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
        'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
        'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl',
        'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
        'Pa','U' ]

class linelist():
  def __init__(self, df, filename=None):
    self.df = df
    self.filename = filename

  @classmethod
  def from_ts_format(self, filename):
    llist = pd.DataFrame()
    with open(filename) as f:
      lines = f.readlines()
      names = ['lam0', 'evi', 'loggf', 'gv', 'upper_g', 'ga', 'gq']
      ncol_max = len(names)

      while len(lines) > 0:
        line = lines.pop(0).strip()
        mass_id, ion, nlines = split_words(line)
        line = lines.pop(0).strip()
        species = line.strip("'").strip()[:6]
        nlines = int(nlines)
        table = ''.join(lines[:nlines])
        df = pd.read_fwf(StringIO(table), header=None)
        del lines[:nlines]
        ncolumns = min(ncol_max, len(df.columns))
        df = df.iloc[:, :ncolumns]
        df.columns = names[:ncolumns]
        if 'gq' in df.columns:
          if df.gq.dtype != float:
            df.drop(columns='gq', inplace=True)
        df['ion'] = species.strip()
        df['iso'] = int(mass_id.strip("'").split('.')[-1])
        llist = pd.concat([llist, df])
    return linelist(llist, filename=filename)
  
  @classmethod
  def from_dispatch_format(self, filename):
    llist = pd.read_csv(filename, sep=r'\s+', quotechar="'")
    llist['lam0'] = vacuum2obs(llist.lam_vac)
    llist = llist.sort_values(by='lam0')
    return linelist(llist, filename=filename)
  
  def write_to_ts_format(self, filename):
    df = self.df
    species = self.df[['ion', 'iso']].drop_duplicates().values
    head = ['lam0', 'evi', 'loggf', 'gv', 'ga', 'gq']
    if hasattr(df, 'label'): head.append('label')
    with open(filename, 'w+') as f:
      for s, iso in species:
        sub = df[(df.ion == s) & (df.iso == iso)]
        nlines = sub.shape[0]
        if len(s.split()) == 2:
          id, ion = s.split()
          if ion=='I': ion = 1
          if ion=='II': ion = 2
          if ion=='III': ion = 3
          if ion=='IV': ion = 4
          if ion=='V': ion = 5
          if ion=='VI': ion = 6
        else:
          id = s
          ion = 1

        if id in element_ids:
          p = element_ids.index(id) + 1
        else:
          p = 0

        f.write(f"'{p}.{iso}' {ion} {nlines}\n")
        f.write(f"'{s}'\n")
        f.write(sub[head].to_string(index=False, header=False) + '\n')