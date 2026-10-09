"""Object views on the tables of a model atom (:mod:`tumult.atom`).

:class:`Rows` is the sequence behind ``atom.levels``, ``atom.bbtrans``, ``atom.bftrans``
and ``atom.colls``; its items are :class:`Row` objects that read the table's columns as
attributes. A bound-free transition spans the cross-section rows of one continuum, so
``atom.bftrans`` groups the rows of ``bft``; the other views have one object per row.
"""
import functools
from collections.abc import Sequence

import pandas as pd

from .const import CC, EE, HH


class Rows(Sequence):
    """The rows of one of an atom's tables as ``kind(rows, k)`` objects, made when
    indexed, or with ``by`` the groups of rows that share the ``by`` columns, in order
    of first appearance. They read the table's columns as numpy arrays, taken when the
    sequence is made: a pandas lookup per attribute would make a table of 10^4
    transitions slow to sort, and an item reflects the table as it was when the view
    was read. The sequence displays as the atom's table and passes other attributes
    (``query``, ``loc``, a column, ...) on to it."""
    def __init__(self, atom, kind, table, by=None):
        df = getattr(atom, table, None)
        self.atom, self.kind, self.table, self.by = atom, kind, table, by
        self.columns = {} if df is None else {c: df[c].to_numpy() for c in df.columns}
        if df is None:
            self.groups = []
        elif by is None:
            self.groups = range(len(df))
        else:
            self.groups = list(df.groupby(by, sort=False).indices.values())

    @functools.cached_property
    def levels(self):
        """The atom's levels, for the ends of its transitions."""
        return self.atom.levels

    def column(self, attr, k):
        """Column ``attr`` of item ``k``: one value per row, and for grouped rows one
        array per group, except a single value for the ``by`` columns."""
        values = self.columns[attr][self.groups[k]]
        return values[0] if self.by is not None and attr in self.by else values

    def frame(self, k):
        """The table's row, or rows, of item ``k``."""
        return getattr(self.atom, self.table).iloc[self.groups[k]]

    def __len__(self):
        return len(self.groups)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self.kind(self, k) for k in range(len(self))[key]]
        k = range(len(self))[key]          # negative indices and IndexError as for a list
        return self.kind(self, k)

    def __getattr__(self, attr):
        table = self.__dict__.get('table')     # absent while a copy is being built
        if table is None:
            raise AttributeError(attr)
        return getattr(getattr(self.atom, table), attr)

    def __repr__(self):
        return getattr(self.atom, self.table).__repr__()

    def _repr_html_(self):
        return getattr(self.atom, self.table)._repr_html_()


class Row(object):
    """Item ``k`` (from 0) of ``rows``, with the table's columns as attributes."""
    def __init__(self, rows, k):
        self.rows = rows
        self.k = k

    def __getattr__(self, attr):
        rows = self.__dict__.get('rows')       # absent while a copy is being built
        if rows is not None and attr in rows.columns:
            return rows.column(attr, self.k)
        raise AttributeError(f'{type(self).__name__} has no attribute {attr!r}')

    def __repr__(self):
        return self.rows.frame(self.k).__repr__()

    def _repr_html_(self):
        frame = self.rows.frame(self.k)
        return frame._repr_html_() if isinstance(frame, pd.DataFrame) else None

    @property
    def atom(self):
        return self.rows.atom


class Level(Row):
    """Energy level ``i`` (from 0) of an atom, row ``i`` of its ``lvl`` table
    (``atom.levels``)."""
    @property
    def i(self):
        return self.k

    @property
    def ev(self):
        return float(self.rows.column('ev', self.k))

    @property
    def ecm(self):
        """The level energy in cm^-1."""
        return self.ev * EE / (CC * HH)

    @property
    def g(self):
        return float(self.rows.column('g', self.k))

    @property
    def ion(self):
        """The ionisation stage, 1 for neutral."""
        return int(self.rows.column('ion', self.k))

    @property
    def label(self):
        """The label the model atom gives the level."""
        return str(self.rows.column('label', self.k))


class LevelPair(Row):
    """A row whose ``i`` and ``j`` columns name a lower and an upper level."""
    @property
    def lower(self):
        """The lower level, a :class:`Level`."""
        return self.rows.levels[int(self.i)]

    @property
    def upper(self):
        """The upper level, a :class:`Level`."""
        return self.rows.levels[int(self.j)]


class BFTransition(LevelPair):
    """Bound-free transition ``k`` (from 0) of an atom (``atom.bftrans``): the rows of
    its ``bft`` table that share ``i`` and, where the table has one, ``j``. ``i`` and
    ``j`` are single values, ``lam_air`` and ``alpha`` the cross section's nodes."""
    @property
    def upper(self):
        """The upper level, a :class:`Level`: the ``j`` column, or without one the
        ground level of the next ionisation stage (``HOP/atom_convert_mod.f90``)."""
        if 'j' in self.rows.columns:
            return super().upper
        levels = self.rows.levels
        stage = self.lower.ion + 1
        return next(lvl for lvl in levels[int(self.i) + 1:] if lvl.ion == stage)


class Collision(LevelPair):
    """Collision record ``k`` (from 0) of an atom, row ``k`` of its ``col`` table
    (``atom.colls``): its ``key``, the levels ``i`` and ``j`` it couples and one column
    per temperature."""
    def __repr__(self):
        # the row mixes text and numbers, and pandas prints such a row's rates of
        # order 1e-9 as 0.0 unless given a float format
        return self.rows.frame(self.k).to_string(float_format='{:.4e}'.format,
                                                 name=True, dtype=True)
