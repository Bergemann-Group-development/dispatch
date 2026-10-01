from io import StringIO
from functools import lru_cache
import logging
from .auxfuncs import split_words
from tumult.const import EK
import os
import re
from tumult.isotopes import parse_ts, format_molec, molec_to_ts, species_elements, element_z
import numpy as np
import pandas as pd
from astropy import units as u
from pathlib import Path

element_ids = [
        'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,
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

ion_dict = {0:'', 1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI', 7:'VII', 8:'VIII', 9:'IX', 10:'X'}
ion_inv = {v: k for k, v in ion_dict.items()}

mol_dict = {101:'H2', 106:'CH', 107:'NH', 108:'OH',
            606:'C2', 707:'N2', 808:'O2', 607:'CN',
            708:'NO', 608:'CO', 10108:'H2O'}

formatters = {'grad':'{:.3e}'.format,
              'gvdw':'{:.8g}'.format,
              'gstk':'{:.3g}'.format,
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
  def from_file(cls, filename, fmt=None):
    if fmt == 'ts':
      return cls.from_ts_format(filename)
    elif fmt == 'tumult':
      return cls.from_tumult_format(filename)
    elif fmt == 'linemake':
      return cls.from_linemake_format(filename)
    elif fmt == 'moog':
      return cls.from_moog_format(filename)
    elif fmt == 'vald':
      return cls.from_vald_format(filename)
    elif fmt == 'abslin':
      return cls.from_abslin_format(filename)
    else:
      raise ValueError(f"Unknown line list format: {fmt}")
    
  @classmethod
  def from_folder(cls, folder, fmt=None, pattern='*'):
    llist = pd.DataFrame()
    for filename in Path(folder).glob(pattern):
      if filename.is_file():
        try:
          llist = pd.concat([llist, cls.from_file(filename, fmt=fmt)])
        except Exception as e:
          print(f"Error reading line list {filename}: {e}")
    return cls(llist)

  @classmethod
  def from_ts_format(cls, filename):
    llist = pd.DataFrame()
    with open(filename) as f:
      lines = f.readlines()
      names = ['lam0', 'evi', 'loggf', 'gvdw', 'upper_g', 'grad', 'gstk']
      levels = ['levlo', 'levup']
      df_dict = {}

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
          # Tokenise directly rather than going through the CSV reader: trailing
          # labels may carry an unbalanced apostrophe from term notation such as
          # "p4f'[5]", which breaks quotechar parsing and silently drops blocks.
          rows = [r.split() for r in lines[:nlines] if r.strip()]

          def isnum(tok):
            try:
              float(tok)
              return True
            except ValueError:
              return False

          # gstk (gamst) is optional, so the designations sit in column 6 or 7
          ncore = 0
          for tok in rows[0][:len(names)]:
            if not isnum(tok):
              break
            ncore += 1

          cols = names[:ncore]
          has_lev = len(rows[0]) >= ncore + 2 and all(
              len(rows[0][ncore + k].strip("'")) == 1 for k in (0, 1))
          if has_lev:
            cols = cols + levels

          data = [r[:len(cols)] for r in rows]
          df = pd.DataFrame(data, columns=cols)
          df[cols[:ncore]] = df[cols[:ncore]].astype(float)
          for lev in levels:
            if lev in df.columns:
              df[lev] = df[lev].str.strip("'").str.strip()

          # decode the Turbospectrum species token into a molecule_mod label:
          # a molecule keeps the sorted, glued mass numbers (each floored to 1
          # so a natural '000' constituent stays present); an atom keeps its
          # bare mass number as an int, as the rest of the code expects
          cons = parse_ts(mass_id)
          if len(cons) > 1:
            iso = format_molec([(z, max(1, a)) for z, a in cons])
          else:
            iso = cons[0][1]

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
        # strength and reach [Angstrom] are the trailing columns linelist_out.txt and
        # hop_prune.x write: the line's peak over continuum extinction, and the half-width
        # out to which its extinction exceeds HOP's tolerance times the continuum
        names = ['lam0', 'evi', 'loggf', 'gvdw', 'grad', 'gstk', 'i', 'j', 'strength', 'reach']
        df_dict = {}
        usecols = lambda col: col is None or col < len(names)

        while len(lines) > 0:
            header = lines.pop(0).strip()
            
            try:
              # species is the first single-quoted field ("'Fe I'"), optionally followed by
              # a label; split on the quote so the ion stage that .split()[0] dropped is kept
              species = header.split("'")[1].strip() if "'" in header else header.split()[0].strip()
              line = lines.pop(0).strip()
              iso_id, nlines = split_words(line)
              iso = int(iso_id.strip("'").split('.')[-1])
              nlines = int(nlines)
              table = ''.join(lines[:nlines])
            except Exception as e:
              print(f"Error parsing line list header in {filename}: {header}")
              raise e

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
  def from_linemake_format(cls, filename):
    llist = pd.read_fwf(filename, widths=[10,10,10,10,20,20], sep=r'\s+', header=None, dtype={1:str})
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
  
  @classmethod
  def from_moog_format(cls, filename):
    names = ['lam0', 'ion', 'evi', 'loggf', 'label']
    df = pd.read_csv(filename, sep=r'\s+', header=None, names=names, dtype={'ion':str})
    mask = df.lam0 > 0
    df = df[mask]
    df['iso'] = 0
    for species in np.unique(df.ion):
      mask = df.ion == species
      p, suf = species.split('.')
      ion = int(suf[0])+1
      iso = int(suf[1:])
      atom = element_ids[int(p)-1]
      df.loc[mask, 'ion'] = atom + ' ' + ion_dict[ion]
      df.loc[mask, 'iso'] = iso

    return cls(df)
  
  @classmethod
  def from_vald_format(cls, filename):
    names = ['ion', 'lam0', 'evi', 'loggf', 'grad', 'gstk', 'gvdw', 'lande', 'label']
    df = pd.read_csv(filename, skiprows=2, header=None, names=names, quotechar="'")
    eid, j = df.ion.str.split(expand=True).values.T
    df['ion'] = [f'{e} {ion_dict[int(i)]}' for e, i in zip(eid, j)]
    df['iso'] = 0
    return cls(df)
  
  @classmethod
  def from_abslin_format(cls, filename):
    df = pd.DataFrame(columns=['ion', 'lam0', 'evi', 'loggf', 'grad', 'gvdw', 'gstk'])

    with open(filename) as f:
      lines = f.readlines()
      lines = [l.strip() for l in lines if len(l.strip()) > 0]
      nlines = int(lines.pop(0))
      while len(lines) > 0:
        # lower level
        idi = lines.pop(0)
        ji, ecmi, gi = [float(x) for x in lines.pop(0).split()]
        # upper level
        idj = lines.pop(0)
        jj, ecmj, gj = [float(x) for x in lines.pop(0).split()]
        # line data
        lam0, lams, lame = [float(x) for x in lines.pop(0).split()]
        fosc, grad, gvdw, gstk = [float(x) for x in lines.pop(0).split()]
        
        evi = (ecmi * u.cm**-1).to(u.eV, equivalencies=u.spectral())
        loggf = np.log10(gi * fosc)

        df.loc[len(df)] = [idi, lam0, evi.value, loggf, grad, gvdw, gstk]

    df.ion = df.ion.str.split().apply(lambda x: x[0].lower().capitalize() + ' ' + x[1])
    df['iso'] = 0
    return cls(df)
  
  def write_to_tumult_format(self, filename, label=None):
    # work with a copy so slice/view issues don't matter
    df = self.copy()
    head = ['lam0', 'evi', 'loggf', 'gvdw', 'grad', 'gstk', 'i', 'j', 'label']
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
    head = ['lam0', 'evi', 'loggf', 'gvdw', 'upper_g', 'grad', 'gstk']
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

        # build the Turbospectrum species token. Atom vs molecule is decided by
        # parsing the species formula against the element table (species_elements
        # falls back to a general chemical-formula parser), so no hardcoded list
        # of molecules is needed -- CaH, SiS, VO, ZrO, ... all resolve on their
        # own. Atoms keep the bare '<Z>.<iso>' form; molecules get the full
        # atomic-number/mass-number encoding.
        elements = species_elements(el_id)
        if len(elements) == 1:
          token = f"{element_z(elements[0])}.{iso}"
        else:
          token = molec_to_ts(el_id, iso)

        f.write(f"'{token}' {ion} {nlines}\n")
        f.write(f"'{s}'\n")
        f.write(sub[head].to_string(index=False, header=False, na_rep='', formatters=formatters) + '\n')


# --- the master line list -------------------------------------------------------
def strip_ion_stage(species):
    """A species label without its ionisation stage: ``"Fe I"`` -> ``"Fe"``.

    A line list names a *species* — a particular ion of an element ("Fe I", "Fe II") or a
    molecule ("CH", "CO"). Dropping the stage leaves the element for an atomic species and
    the molecule unchanged for a molecular one, which carries no stage: ``"CH"`` -> ``"CH"``.
    So the result is an element symbol only when the input was an atomic species; for the
    constituent elements of a molecule see :func:`tumult.isotopes.species_elements`.

    Whitespace is normalised first, so ``"Fe  II"`` and ``"Fe II"`` answer alike.
    """
    parts = str(species).split()
    return parts[0] if parts else ""


def strip_ion_stages(species):
    """:func:`strip_ion_stage` over a line list's ``ion`` column, as a pandas Series."""
    return species.astype(str).str.split().str[0]


MASTER_PATH = os.environ.get(
    "TUMULT_MASTER_LINELIST",
    str(Path(__file__).resolve().parent.parent / "input" / "linelists" / "master_line.txt"),
)

_DUMMY_T = 5000.0   # excitation temperature of the strength proxy below
_cache = {}  # path -> parsed DataFrame (process-wide, read-only, shared across sessions)


def species_header(path):
    """The first quoted species label in a Tumult line-list file."""
    with Path(path).open(errors="ignore") as stream:
        header = stream.readline().strip()
    if "'" in header:
        return header.split("'", 2)[1].strip()
    return header.split()[0] if header.split() else ""


def is_molecule(species):
    """Whether a species label has more than one constituent; unknown labels are kept."""
    try:
        return len(species_elements(species)) > 1
    except (IndexError, ValueError):
        return False


# composition_patch_mod allocates this many linelist_file slots
MAX_LINELIST_FILES = 128


def linelist_files(path):
    """The line-list files a file-or-folder setting names, as resolved paths: the file
    itself, or a folder's direct files in sorted order (dotfiles left out). A folder is
    expanded here rather than handed to dispatch as ``linelist_folder``: its listing
    there goes through a shell temp file in the working directory, which the concurrent
    runs of a synthesis clobber for each other, and the order is then the shell's."""
    path = Path(path).expanduser().resolve()
    if not path.is_dir():
        return (path,)
    files = tuple(sorted(item for item in path.iterdir()
                         if item.is_file() and not item.name.startswith('.')))
    if len(files) > MAX_LINELIST_FILES:
        raise ValueError(f"{path} holds {len(files)} line-list files; dispatch reads at "
                         f"most {MAX_LINELIST_FILES}")
    return files


@lru_cache(maxsize=8)
def atomic_files(folder):
    """Files in ``folder`` whose first species is atomic.

    A conventional atomic label has one constituent (``Fe I`` -> Fe), whereas molecular
    labels have several (``TiO`` -> Ti + O). Files with an unfamiliar header stay in the
    selection: omitting an unclassified input would silently lose annotations.
    """
    files = []
    for path in sorted(item for item in Path(folder).iterdir()
                       if item.is_file() and not item.name.startswith(".")):
        try:
            molecular = is_molecule(species_header(path))
        except OSError:
            logging.warning("Could not classify line-list file %s; keeping it", path)
            molecular = False
        if not molecular:
            files.append(path)
    return tuple(files)


HLINEDATA_NAME = "Hlinedata"   # HOP's dedicated hydrogen file (composition_mod add_h_block)


def _read_hlinedata(path):
    """The dedicated hydrogen file in the frame layout of :meth:`LineList.from_tumult_format`.

    Its columns are ``lam_air, i, j, lower_ev, upper_ev, loggf, label, grad``; HOP selects
    this layout by the file name, so the generic reader must not see it.
    """
    with Path(path).open() as stream:
        stream.readline()                      # the descriptive header
        nlines = int(split_words(stream.readline())[1])
        table = "".join(stream.readline() for _ in range(nlines))
    df = pd.read_csv(StringIO(table), header=None, sep=r"\s+", quotechar="'",
                     names=["lam0", "i", "j", "evi", "eup", "loggf", "label", "grad"])
    df["ion"], df["iso"] = "H I", 0
    return df[["ion", "iso", "lam0", "evi", "loggf", "i", "j", "grad"]]


def _read_list(path):
    path = Path(path)
    if path.name == HLINEDATA_NAME:
        return _read_hlinedata(path)
    return LineList.from_tumult_format(path).reset_index()


def load_master(path=MASTER_PATH):
    """Parse a line list once per path, ensuring it carries a ``strength`` column.

    ``path`` is a file or a folder of per-species files, the same input the synthesis
    takes. A folder is read as its atomic files plus :data:`HLINEDATA_NAME`: the whole
    frame is held in memory for searching, and the molecular lists (hundreds of MB for
    the default folder) belong to the synthesis, which crops them to its windows.

    A ``linelist_out.txt`` written by a run already has one — HOP's line-centre over
    continuum extinction, which compares across species — and it is kept as it is. Only a
    plain list without that column gets the Boltzmann proxy below, which is an ordering of
    last resort: it has no abundance, no ionization fraction and no partition function in
    it, so it is meaningful within one ion and misleading between two. ``attrs`` records
    which of the two happened, and :mod:`linestrength` labels the table accordingly.
    """
    path = str(path)
    if path not in _cache:
        source = Path(path)
        files = atomic_files(source) if source.is_dir() else (source,)
        if not files:
            raise ValueError(f"No line lists in {source}")
        m = pd.concat([_read_list(f) for f in files], ignore_index=True)
        proxy = "strength" not in m.columns or not m.strength.notna().any()
        if proxy:
            m["strength"] = m.lam0 * np.exp(-m.evi * EK / _DUMMY_T) * 10.0 ** m.loggf
        m.attrs["strength_is_proxy"] = proxy
        _cache[path] = m
    return _cache[path]


def forget(path):
    """Drop the parsed frame of ``path`` (a file or folder) and the folder
    classification, so the next :func:`load_master` reads the disk again; for a folder
    whose files were added or removed."""
    _cache.pop(str(path), None)
    atomic_files.cache_clear()


# --- searching a line list ---------------------------------------------------------
def search_lines(lines, query):
    """The rows of ``lines`` matching a free-text ``query`` such as ``"O I 77"``.

    The query is split on whitespace. Purely numeric words (``77``, ``7771.9``) match the
    wavelength ``lam0`` in Angstrom, everything else names the species: those words, joined
    by single spaces, are matched word by word against the ``ion`` label, case-insensitively,
    the leading words exactly and the last one as a prefix. So ``"O"`` matches every species
    starting with O (O I, O II, OH, Os I), ``"O I"`` only O I and ``"Fe"`` both Fe I and Fe
    II. The ionisation stage is matched exactly once a further word follows (``"O I "``).

    A numeric word matches as the digits are typed: ``77`` is the leading digits of any
    wavelength, so it finds 77xx and 77xxx Å alike, and ``777`` narrows to 777x. The
    wavelength may be read off rounded as well as truncated (``7772`` finds 7771.959).
    A decimal point pins the magnitude. Several numeric words match any of them.

    Returns a DataFrame (empty for an empty query); ``lines`` must carry ``ion`` and
    ``lam0`` columns.
    """
    words = str(query).split()
    if not words:
        return lines.iloc[0:0]
    numeric = [w for w in words if _NUMERIC.fullmatch(w)]
    species = [w for w in words if not _NUMERIC.fullmatch(w)]
    # the species is still being typed only while it is the last thing in the query
    complete = words[-1] not in species or str(query)[-1].isspace()
    mask = np.ones(len(lines), dtype=bool)
    if species:
        mask &= _species_mask(lines.ion, species, complete)
    if numeric:
        lam = lines.lam0.values.astype(float)
        hit = np.zeros(len(lines), dtype=bool)
        for word in numeric:
            hit |= _wavelength_mask(lam, word)
        mask &= hit
    return lines[mask]


_NUMERIC = re.compile(r"\d+\.?\d*|\.\d+")


def _species_mask(ion, words, complete):
    """Rows whose ion label matches the species ``words`` (see :func:`search_lines`);
    with ``complete`` the last word is matched exactly rather than as a prefix.

    Matched over the distinct labels rather than the column: a master list has a few
    hundred species over a few hundred thousand rows, and this runs on every keystroke."""
    words = [w.lower() for w in words]

    def matches(label):
        parts = str(label).lower().split()
        if len(parts) < len(words):
            return False
        head, last = words[:-1], words[-1]
        return (parts[:len(head)] == head
                and (parts[len(head)] == last if complete
                     else parts[len(head)].startswith(last)))

    labels = [lab for lab in ion.unique() if matches(lab)]
    return ion.isin(labels).values


def _wavelength_mask(lam, word):
    """Rows whose ``lam`` reads as ``word`` when written out, truncated or rounded, to the
    precision the word is typed at; an integer word may be the leading digits of a longer
    wavelength. Wavelengths in Å, from 1 Å to below 1e7 Å."""
    value = float(word)
    if "." in word:
        step = 10.0 ** -(len(word) - word.index(".") - 1)
        shifts = (1.0,)
    else:
        step = 1.0
        shifts = tuple(10.0 ** k for k in range(7 - len(word)))
    mask = np.zeros(len(lam), dtype=bool)
    for shift in shifts:
        lo = value * shift - step * shift / 2 if shift == 1.0 else value * shift
        # rounding at the typed precision reaches half a step below the word's value; a
        # shifted (longer) wavelength is only ever read by its leading digits
        mask |= (lam >= lo) & (lam < (value + step) * shift)
    return mask
