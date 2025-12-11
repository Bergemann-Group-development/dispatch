
from requests import head
from tumult.auxfuncs import *
import numpy as np
import pandas as pd
import textwrap

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

ion_dict = {0:'', 1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI'}
ion_inv = {v: k for k, v in ion_dict.items()}

mol_dict = {101:'H2', 106:'CH', 107:'NH', 108:'OH',
            606:'C2', 707:'N2', 808:'O2', 607:'CN',
            708:'NO', 608:'CO', 10108:'H2O'}

formatters = {'ga':'{:.3e}'.format,
              'gv':'{:.8g}'.format,
              'gq':'{:.3g}'.format,
              'ev':'{:.6f}'.format,
              'evi':'{:.3f}'.format,
              'lam_air':'{:.3f}'.format,
              'lam0':'{:.3f}'.format,
              'alpha':'{:.3e}'.format,
              'label':"'{:}'".format,
              'loggf':'{:.3f}'.format}

class LineList(pd.DataFrame):
  # Required for pandas to return instances of LineList on slicing, etc.
  @property
  def _constructor(self):
      return LineList

  @classmethod
  def from_ts_format(cls, filename):
    llist = pd.DataFrame()
    with open(filename) as f:
      lines = f.readlines()
      names = ['lam0', 'evi', 'loggf', 'gv', 'upper_g', 'ga', 'gq']
      df_dict = {}
      usecols = lambda col: col is None or col < len(names)

      while len(lines) > 0:
        try:
          line = lines.pop(0).strip()
          mass_id, ion, nlines = split_words(line)
        except Exception as e:
          print(f"Error parsing line list header in {filename}: {line}")
          raise e
        
        line = lines.pop(0).strip()
        label = line.strip("'").strip()
        label = ' '.join(label.split())
        nlines = int(nlines)
        table = ''.join(lines[:nlines])

        try:
          try:
            df = pd.read_csv(StringIO(table), header=None, sep=r'\s+', quotechar="'", usecols=usecols)
          except Exception as e:
            print(f"Retrying with fixed-width format for {label} due to error: {e}")
            df = pd.read_fwf(StringIO(table), header=None, infer_nrows=nlines).iloc[:, :len(names)]
          
          df.columns = names[:df.shape[1]]

          if 'gq' in df.columns:
            if df.gq.dtype != float:
              df.drop(columns='gq', inplace=True)

          npro = mass_id.strip("'").split('.')[0]
          iso = mass_id.strip("'").split('.')[-1]
          
          if int(npro) > 100:
            iso = textwrap.wrap(iso, 3)
            iso = sorted([max(1, int(i)) for i in iso])
            iso = ''.join([str(i) for i in iso])
          else:
            iso = int(iso)

          df_dict[(label, iso)] = df
        except Exception as e:
          print(f"Error reading line list {filename} for {label}: {e}")

        del lines[:nlines]
      
      llist = pd.concat(df_dict, names=['ion', 'iso', None])
      llist = llist.sort_index()

    return cls(llist)
  
  @classmethod
  def from_tumult_format(cls, filename):
    llist = pd.DataFrame()
    with open(filename) as f:
        lines = f.readlines()
        names = ['lam0', 'evi', 'loggf', 'gv', 'ga', 'gq', 'i', 'j']
        df_dict = {}
        usecols = lambda col: col is None or col < len(names)

        while len(lines) > 0:
            species = lines.pop(0).strip()
            line = lines.pop(0).strip()
            iso_id, nlines = split_words(line)
            iso = int(iso_id.strip("'").split('.')[-1])
            nlines = int(nlines)
            table = ''.join(lines[:nlines])

            try:
              df = pd.read_csv(StringIO(table), header=None, sep=r'\s+', quotechar="'", usecols=usecols)
            except Exception as e:
              print(f"Retrying with fixed-width format for {species} due to error: {e}")
              df = pd.read_fwf(StringIO(table), header=None, infer_nrows=nlines).iloc[:, :len(names)]
              
            df.columns = names[:df.shape[1]]
            del lines[:nlines]

            ion = species.strip()
            df_dict[(ion, iso)] = df
        
        llist = pd.concat(df_dict, names=['ion', 'iso', None])
        llist = llist.sort_index()

    return cls(llist)
  
  @classmethod
  def from_dispatch_format(cls, filename):
    llist = pd.read_csv(filename, sep=r'\s+', quotechar="'")
    llist = llist.rename(columns={'lam_air':'lam0'})
    llist = llist.sort_values(by='lam0')
    return cls(llist)
  
  @classmethod
  def from_linemake(cls, file):
    llist = pd.read_fwf(file, widths=[10,10,10,10,20,20], sep=r'\s+', header=None, dtype={1:str})
    llist.rename(columns={0: 'lam0', 1: 'id', 2: 'evi', 3: 'loggf', llist.shape[1]-1: 'label'}, inplace=True)
    llist = llist[['lam0', 'id', 'evi', 'loggf', 'label']]
    llist[['np', 'suf']] = llist.id.str.split('.', expand=True)
    atom = llist.np.astype(int) < len(element_ids)
    mols = llist.np.str.contains('0') & ~atom
    llist.loc[mols, 'id'] = np.vectorize(mol_dict.get)(llist.np[mols].astype(int).values)

    llist['iso'] = llist.suf.str[1:]
    llist.loc[llist.iso=='', 'iso'] = '0'
    llist.loc[atom, 'ion'] = llist.suf[atom].str[:1].astype(int) + 1
    llist.loc[atom, 'ion_id'] = np.vectorize(ion_dict.get)(llist.ion[atom].values)
    llist.loc[atom, 'id'] = np.array(element_ids)[llist.np.astype(int)[atom]-1]
    llist.loc[atom, 'id'] = llist.loc[atom, ['id', 'ion_id']].agg(' '.join, axis=1)
    llist.drop(columns=['ion_id', 'ion', 'np', 'suf'], inplace=True)
    llist.rename(columns={'id': 'ion'}, inplace=True)
    llist.loc[llist.label.isna(), 'label'] = ''

    mask = mols | atom
    llist = llist[mask]
    return cls(llist)
  
  def write_to_tumult_format(self, filename, label=None):
    # work with a copy so slice/view issues don't matter
    df = self.copy()
    head = ['lam0', 'evi', 'loggf', 'gv', 'ga', 'gq', 'i', 'j', 'label']
    present = [col for col in head if col in df.columns]
    head = head[:head.index(present[-1]) + 1]
    df = df.reindex(columns=df.columns.union(head), fill_value=0)
    df = df.fillna(0)

    with open(filename, 'w+') as f:
      for (s, iso), sub in df.groupby(['ion', 'iso']):
        nlines = len(sub)
        sub = sub[head].sort_values(by='lam0')

        if f.tell() > 0: f.write('\n')
        f.write(f"'{s}'")
        if label: f.write(f" '{label}'")
        f.write("\n")
        f.write(f"'{iso}' {nlines}\n")
        f.write(sub.to_string(index=False, header=False, na_rep='', formatters=formatters))
  
  def write_to_ts_format(self, filename):
    # work with a copy so slice/view issues don't matter
    df = self.copy()

    # ensure all required columns exist (fill missing with 0)
    head = ['lam0', 'evi', 'loggf', 'gv', 'upper_g', 'ga', 'gq']
    df = df.reindex(columns=df.columns.union(head), fill_value=0)
    df['upper_g'] = np.maximum(df.upper_g, 1.0)  # avoid zero upper_g

    if 'label' in df.columns:
      head.append('label')

    with open(filename, 'w+') as f:
      for (s, iso), sub in df.groupby(['ion', 'iso']):
        nlines = len(sub)

        # parse species string into element + ion
        if len(s.split()) == 2:
          el_id, ion = s.split()
          ion = ion_inv.get(ion, 1)
        else:
          el_id, ion = s, 1

        # map element id to integer code
        p = element_ids.index(el_id) + 1 if el_id in element_ids else 0

        f.write(f"'{p}.{iso}' {ion} {nlines}\n")
        f.write(f"'{s}'\n")
        f.write(sub[head].to_string(index=False, header=False, na_rep='', formatters=formatters) + '\n')