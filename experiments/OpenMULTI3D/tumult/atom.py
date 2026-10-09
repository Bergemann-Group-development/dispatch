"""Model atoms in the table format, the one format dispatch.x reads.

A run converts a classic MULTI atom (a file named ``atom.*``) to the table format with
``HOP/atom_convert_mod.f90`` and keeps the conversion in its ``atoms/`` folder;
:meth:`Atom.from_old` applies the same conversion through ``atom_convert.x``. The tables
are pandas DataFrames: levels (``lvl``), bound-bound transitions (``bbt``), bound-free
cross-sections (``bft``) and collision records (``col``). :attr:`Atom.levels`,
:attr:`Atom.bbtrans`, :attr:`Atom.bftrans` and :attr:`Atom.colls` give objects on top of
them (:mod:`tumult.atom_rows`).
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .auxfuncs import get_qmax
from .const import CC, UU
from .atom_rows import BFTransition, Collision, Level, Rows
from .trans import Transition

CONVERTER = Path(__file__).resolve().parents[1] / 'HOP' / 'standalone' / 'build' / 'atom_convert.x'
MOLECULES = ['H2', 'C2', 'N2', 'O2', 'CH', 'OH', 'NH', 'CN', 'CO', 'NO', 'H2O']

formatters = {'grad':'{:.3e}'.format,
                            'gvdw':'{:.8g}'.format,
                            'gstk':'{:.3g}'.format,
                            'ev':'{:.6f}'.format,
                            'evi':'{:.6f}'.format,
                            'lam_air':'{:.3f}'.format,
                            'alpha':'{:.3e}'.format,
                            'label':"'{:}'".format,
                            'loggf':'{:.3f}'.format}

def read(text, dtype):
    return np.array(next(text).strip().split()).astype(dtype)


class Atom(object):
    """A model atom in the table format; :meth:`read` reads either format from file."""
    def __init__(self):
        pass

    @property
    def levels(self):
        """The energy levels as :class:`Level` objects, indexed from 0."""
        return Rows(self, Level, 'lvl')

    @property
    def bbtrans(self):
        """The bound-bound transitions as :class:`tumult.trans.Transition` objects, in
        table order."""
        return Rows(self, Transition, 'bbt')

    @property
    def bftrans(self):
        """The bound-free transitions as :class:`tumult.atom_rows.BFTransition` objects,
        one per continuum, in table order."""
        by = ['i', 'j'] if 'j' in getattr(self, 'bft', ()) else ['i']
        return Rows(self, BFTransition, 'bft', by=by)

    @property
    def colls(self):
        """The collision records as :class:`tumult.atom_rows.Collision` objects, in table
        order."""
        return Rows(self, Collision, 'col')

    @property
    def ground_states(self):
        """The index of the first level of each ionisation stage."""
        return np.where(np.diff(self.lvl.ion.to_numpy(), prepend=0) > 0)[0]

    @property
    def ncont(self):
        """The number of bound-free transitions (continua)."""
        return len(self.bftrans)

    def plot_grotrian(self, **kwargs):
        plot_grotrian(self.lvl, self.bbt, **kwargs)

    def post_init(self):
        self.molec = self.id in MOLECULES
        if not self.molec: self.id = self.id.capitalize()
        self.awgt = self.ai * UU

        self.lvl_tab = self.lvl.columns
        if hasattr(self, 'bft'): self.bft_tab = self.bft.columns
        if hasattr(self, 'col'): self.col_tab = self.col.columns
        if hasattr(self, 'bbt'): self.bbt_tab = self.bbt.columns

    @classmethod
    def read(cls, atom_file):
        """Read a model atom: a file named ``atom.*`` is a MULTI atom and is converted as
        a run converts it (:meth:`from_old`), any other is in the table format."""
        if Path(atom_file).name.startswith('atom.'):
            return cls.from_old(atom_file)
        return cls.read_table(atom_file)

    @classmethod
    def read_table(cls, atom_file):
        """Read a model atom in the table format."""
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
        if self.nline > 0:
            self.bbt = pd.read_csv(nrows=self.nline, **pdkwargs)
            self.bbt = self.bbt.rename(columns={'lam_air':'lam0'})
            pdkwargs['header'] = 3+self.nk+1+self.nline+1
        if nnu_bf > 0:
            self.bft = pd.read_csv(nrows=nnu_bf, **pdkwargs)
            pdkwargs['header'] = 3+self.nk+1+self.nline+1+nnu_bf+1
        if self.ncol > 0:
            self.col = pd.read_csv(nrows=self.ncol, **pdkwargs)
        self.post_init()
        return self

    @classmethod
    def from_old(cls, atom_file):
        """Convert a classic MULTI atom with ``atom_convert.x``, the conversion dispatch.x
        applies to the MULTI atoms of a run (``HOP/atom_convert_mod.f90``)."""
        if not CONVERTER.is_file():
            raise FileNotFoundError(f'atom_convert.x not found at {CONVERTER}; build it with '
                                    f'"make -C {CONVERTER.parents[1]}"')
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'converted.txt'
            proc = subprocess.run([str(CONVERTER), str(atom_file), str(out)],
                                  capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f'atom_convert.x failed on {atom_file}: '
                                   f'{(proc.stderr or proc.stdout).strip()}')
            return cls.read_table(out)

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

    def get_title(self, abnd=None):
            if abnd is None: abnd = self.abnd
            return 'A({}) = {:.2f}'.format(self.id.strip(), abnd)

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
        """The wavelength range of each transition (``lb``, ``lr`` [AA]) from its profile
        extent, and its indices on the run's grid ``lam`` (``ib``, ``i0``, ``ir``)."""
        self.qnorm = qnorm * 1e5 / CC
        mask = self.bbt.nq > 0
        qmax = get_qmax(self.bbt.qmax[mask], self.bbt.q0[mask]) * self.qnorm
        self.bbt.loc[mask, 'lb'] = self.bbt.loc[mask, 'lam0'] * (1 - qmax)
        self.bbt.loc[mask, 'lr'] = self.bbt.loc[mask, 'lam0'] * (1 + qmax)
        self.bbt.loc[~mask, 'lb'] = self.bbt.loc[~mask, 'lam0'] - self.bbt.loc[~mask, 'qmax']
        self.bbt.loc[~mask, 'lr'] = self.bbt.loc[~mask, 'lam0'] + self.bbt.loc[~mask, 'qmax']
        self.bbt['i0'] = np.searchsorted(lam, self.bbt.lam0)
        self.bbt['ib'] = np.searchsorted(lam, self.bbt.lb)
        self.bbt['ir'] = np.searchsorted(lam, self.bbt.lr)

    def load_non_overlapping(self, nnu, nq_fac):
        """The indices of each transition when the run saves the transitions' own
        frequency grids one after another at the end of its spectrum, in table order
        (``overlapping_bb_in_output`` off): ``nq_fac (nq-1) + 1`` points each, odd nq."""
        nq = self.bbt.nq.to_numpy()
        nq = nq + (1-nq%2)
        nq = nq_fac * (nq-1) + 1
        ir = np.ones_like(nq)*nnu
        ir[1:] = nnu - np.cumsum(nq[:-1])
        ib = ir - nq
        self.bbt['ib'] = ib
        self.bbt['ir'] = ir
        self.bbt['i0'] = (ib + ir) // 2

    def to_linelist(self, file, bbt=None):
        if bbt is None: bbt = self.bbt
        f = open(file, 'w+')
        ions = np.unique(self.lvl.ion).astype(int)
        for ion in ions:
            ll = pd.DataFrame()
            mask = np.array(self.lvl.ion[bbt.i] == ion)
            ev0 = self.lvl.ev[self.lvl.ion == ion].values[0]
            if any(mask):
                f.write(f"'  {np.round(self.ai).astype(int):<16.3f}' {ion:d} {np.sum(mask):d}\n")
                f.write(f"'{self.id:3} {'I'*ion:<3}'\n")
                bb = bbt[mask]
                ll['lam0'] = bb['lam0']
                ll['evi'] = bb['evi'] - ev0
                ll['loggf'] = bb['loggf']
                ll['gvdw'] = bb['gvdw']
                ll['g_j'] = self.lvl.g[bb.j].values
                ll['grad'] = bb['grad']
                ll['gstk'] = bb['gstk']
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


def read_atom(atom_file):
    """Read a model atom in either format (:meth:`Atom.read`)."""
    return Atom.read(atom_file)


def convert_label_to_x_axis(label: str):
    # ok, this is a bit tricky
    # the label is in the form "x4D1*"
    # zyx or abc -> parity. zyx is odd, abc is even
    # 4 is the s-value
    # D is the L-value
    # 1 is the J-value
    # * is the parity odd/even. with * it is odd, without it is even

        label_og = label

        # sometimes however, the label is something like "l4D" or just "l4D*"
        # so we need to handle that. we dont care about J-value actually.
        # we need to return s-value, L-value and * if it is there
        # if length < 3, return None
        if len(label) < 3:
                print(f"Label {label_og} is too short")
                return None

            # first lets cut down to last 4 characters
        label = label[-4:]
        # now lets see the last character
        if label[-1] == "*":
                parity_even = False
                # pop last character
                label = label[:-1]
        else:
                parity_even = True
            # now lets see the last character. if a number, pop it
        if label[-1].isnumeric():
                J_value = int(label[-1])
                label = label[:-1]
            # ok now the last 2 values should hopefully be s-value and L-value
        L_value = label[-1]
        s_value = label[-2]

        # if s_value is not a number, then return None
        if not s_value.isnumeric():
                print(f"Label {label_og} has non-numeric s-value")
                return None
        s_value = int(s_value)

        # if L_value is not a letter, then return None
        if not L_value.isalpha():
                print(f"Label {label_og} has non-letter L-value")
                return None

            # now return the values
        value_to_return = f"{s_value}{L_value}"
        if not parity_even:
                value_to_return += "*"
        return value_to_return

def plot_grotrian(lvl, bbt, ion=None, lw=0.8, level_color='k', level_lw=1, color='lightgrey', noff=2, **kwargs):
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        icont = np.where(np.diff(lvl.ion, prepend=-1)>0)[0]

        if ion is not None:
                lvl = lvl.iloc[:icont[ion]+1]

        ev = lvl.ev.values
        olabels = lvl.label.str.split(expand=True).iloc[:, -1].values
        labels = [convert_label_to_x_axis(label) for label in olabels]
        mask = np.array([label is None for label in labels])
        labels = np.array([label.replace('*', r'$^o$') if label is not None else None for label in labels])

        num_labels = [l
                        .replace('S', '1')
                        .replace('P', '2')
                        .replace('D', '3')
                        .replace('F', '4')
                        if l is not None else '0'
                        for l in labels]
    
        labels[mask] = olabels[mask]
        ulabels, idx, xlabel = np.unique(num_labels, return_index=True, return_inverse=True)

        lower = bbt.i
        upper = bbt.j
        lvlmin = lvl.index[0]
        lvlmax = lvl.index[-1]
        mask = (lower >= lvlmin) & (upper <= lvlmax)
        lower = lower[mask] - lvlmin
        upper = upper[mask] - lvlmin

        x = [xlabel[lower], xlabel[upper]]
        y = [ev[lower], ev[upper]]

        xlvl = [xlabel-0.2, xlabel+0.2]
        ylvl = [ev, ev]

        ax = plt.gca()
        ax.plot(x, y, color=color, lw=lw, **kwargs)
        ax.plot(xlvl, ylvl, color=level_color, lw=level_lw)

        if ion is not None:
                iconts = icont[ion:ion+1]
        else:
                iconts = icont
        
        for i in iconts:
                ax.axhline(lvl.ev[i], color='lightgrey', lw=1, ls='--', zorder=-10)

        ax.set_xticks(np.arange(len(ulabels)))
        ax.set_xticklabels(labels[idx])
        ax.tick_params(axis='x', which='minor', length=0)
        ax.invert_xaxis()
        ax.set_ylabel('Level energy [eV]')
        #title = 'Grotrian Diagram of {:} '.format(self.id)
        #if ion is not None: title = title + 'I'*ion
        #ax.set_title(title)

        plabels = ax.get_xticklabels()
        # Loop over the labels and adjust the y position for every second label
        for i, label in enumerate(plabels):
                off = -0.05*(i%noff)
                label.set_y(off)  # adjust the position

