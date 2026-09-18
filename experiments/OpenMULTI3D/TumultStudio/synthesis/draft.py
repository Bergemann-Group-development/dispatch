"""Calculation drafts: local edits become last-used preferences only on submission."""
from copy import deepcopy

from ..workspace import store


class CalculationDraft:
    def __init__(self, forms=(), setup=None, atoms=None, owner=None, attributes=()):
        self.forms = forms  # (preference section, Form)
        self.setup = setup
        self.atoms = atoms
        self.owner = owner
        self.attributes = attributes
        self.baseline = None

    def begin(self):
        if self.baseline is not None:
            return
        self.baseline = ([deepcopy(form.values()) for _, form in self.forms],
                         {key: deepcopy(getattr(self.owner, key)) for key in self.attributes},
                         list(self.atoms.selected) if self.atoms else None)

    def commit(self):
        if self.baseline is None:
            return
        if self.setup is not None:
            self.setup.commit()
        with store.transaction() as data:
            for section, form in self.forms:
                data['preferences'].setdefault(section, {}).update(form.values())
            if self.atoms is not None:
                data['preferences'].setdefault(self.atoms.section, {})[self.atoms.key] = list(self.atoms.selected)
        self.baseline = None

    def discard(self):
        if self.baseline is None:
            return
        values, attrs, atoms = self.baseline
        for (_, form), saved in zip(self.forms, values):
            form.set_values(**saved)
        for key, value in attrs.items():
            setattr(self.owner, key, value)
        if self.setup is not None:
            self.setup.discard()
        if self.atoms is not None:
            self.atoms.selected = atoms
            self.atoms.render()
        self.baseline = None
