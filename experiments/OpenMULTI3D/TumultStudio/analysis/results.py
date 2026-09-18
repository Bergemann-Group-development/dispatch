"""The analysis results of one document: the collections, their identities and the GUI
half of staleness.

Each of the three analysis tools owns one :class:`ResultCollection` — an ordered list of
result objects (:class:`tumult.abundance.AbundanceResult`,
:class:`tumult.measurement.BalanceResult`) and the uid of the one in view. A result is
identified by its ``uid`` alone, everywhere: the sidebar row, the workflow record that
continues it, the Edit menu, the removal/Undo path and the saved workspace. Names are for
people and are made unique by the document; nothing is ever matched by one.

Saved workspaces are read back here too, as plain functions rather than through the views:
:func:`deserialise` builds a section's result objects and :func:`relink` clears the
references that name objects the load does not restore, so ``workspace.session`` can
interpret a slot without constructing a single Bokeh view (see
:func:`tumult.abundance.AbundanceResult.from_dict` for the object side of it).

Staleness has two halves and this module owns the GUI one. Whether a solution still matches
its *scientific* inputs is the result type's own question (``res.stale_reasons(digest)``);
what this adds is the lookup only the document can do — which loaded spectrum a result's
``ref_uid`` names, that spectrum's digest as it stands now, and how the preparation it was
measured under differs from the preparation that spectrum currently carries.
"""
from __future__ import annotations

import uuid
from copy import deepcopy

from tumult import abundance as ab
from tumult import measurement as ms
from tumult.measured import target_digest

DIGEST_CACHE = 32            # spectrum objects whose digest is remembered at a time


class ResultCollection:
    """One tool's results and the uid of the selected one.

    The list order is the sidebar order and the insertion order of the document; it is what
    a saved workspace records and restores. Nothing here draws or notifies: the view that
    owns the collection decides what a change repaints.
    """

    def __init__(self):
        self.results = []
        self.selected = None

    def __iter__(self):
        return iter(self.results)

    def __len__(self):
        return len(self.results)

    def __getitem__(self, index):
        return self.results[index]

    def find(self, uid):
        return next((r for r in self.results if r.uid == uid), None)

    def selected_result(self):
        return self.find(self.selected)

    def index(self, res):
        return self.results.index(res)

    def append(self, res):
        """List ``res`` and select it (it is what the document just produced)."""
        self.results.append(res)
        self.selected = res.uid

    def insert(self, res, index):
        """Put ``res`` back at ``index`` (the removal Undo); a duplicate uid is a bug."""
        if self.find(res.uid) is not None:
            raise ValueError(f"Duplicate result identity: {res.uid}")
        self.results.insert(min(index, len(self.results)), res)

    def remove(self, res):
        """Drop ``res``; when it was selected, the last remaining result takes its place.
        Returns whether it was the selected one."""
        was_selected = self.selected == res.uid
        self.results.remove(res)
        if was_selected:
            self.selected = self.results[-1].uid if self.results else None
        return was_selected

    def duplicate(self, res):
        """A copy of ``res`` with an identity of its own, not yet listed: the original's
        workflow keeps reaching the original only, and curating either leaves the other as
        it is."""
        copy = deepcopy(res)
        copy.uid = uuid.uuid4().hex[:8]
        return copy

    def replace(self, results, selected):
        """Take over a restored list; ``selected`` must be one of them or ``None``."""
        self.results = list(results)
        uids = {r.uid for r in self.results}
        self.selected = selected if selected in uids else None

    def snapshot(self):
        """``{"results": [...], "selected": uid}`` for the session slot."""
        return {"results": [r.to_dict() for r in self.results], "selected": self.selected}


class ResultReferences:
    """What a result's recorded inputs mean in *this* document.

    ``spec_of(res)`` resolves a result's reference observation to the loaded spectrum object
    (``None`` once it is gone) and the canvas states the preparation that spectrum currently
    carries. :meth:`reasons` hands the result type the digest of what it would now be solved
    from and adds the preparation differences the result type cannot see.
    """

    def __init__(self, canvas, spec_of):
        self.canvas = canvas
        self.spec_of = spec_of
        self._digests = {}               # id(spectrum) -> (spectrum, target_digest)

    def digest(self, spec):
        """``target_digest(spec)``, remembered per spectrum object: the canvas replaces a
        target's spectrum object when its continuum or filter changes."""
        cached = self._digests.get(id(spec))
        if cached is None or cached[0] is not spec:
            if len(self._digests) > DIGEST_CACHE:
                self._digests.clear()
            cached = self._digests[id(spec)] = (spec, target_digest(spec))
        return cached[1]

    def preparation_changes(self, res):
        """What differs between the preparation ``res`` recorded and its target's current
        one, in words; empty when either is unknown. Only what the measurement reads counts:
        the convolution stays in the record as provenance but is not applied to the target
        the lines are measured on."""
        recorded = res.preparation or {}
        obj = self.canvas.spectra.find(res.ref_uid) if res.ref_uid is not None else None
        if not recorded or obj is None:
            return []
        now = self.canvas.preparation(obj)
        changes = []
        if bool(recorded.get("filtered")) != now["filtered"]:
            changes.append("filtering was switched " + ("on" if now["filtered"] else "off"))
        if recorded.get("normalisation", now["normalisation"]) != now["normalisation"]:
            changes.append(f"the target now uses {now['normalisation']}")
        return changes

    def reasons(self, res):
        """``{key: text}`` saying why the solved ``res`` no longer matches its inputs."""
        spec = self.spec_of(res) if res.ref_uid is not None else None
        reasons = res.stale_reasons(None if spec is None else self.digest(spec))
        changes = self.preparation_changes(res)
        if changes:
            reasons["target"] = "; ".join(changes) + " after this was solved"
        elif "target" in reasons:
            reasons["target"] = ("the target's continuum or filtered flux changed after this "
                                 "was solved")
        return reasons


# --- saved workspaces: the result sections, read back without a view --------------------

def _relink_reference(results, spectrum_uids, label):
    """Clear the reference observation of every result whose spectrum a load does not
    restore; one note per cleared reference."""
    notes = []
    for res in results:
        if res.ref_uid is not None and res.ref_uid not in spectrum_uids:
            notes.append(f"{label} result {res.name or res.uid}: reference spectrum "
                         f"{res.ref_name or res.ref_uid} is not loaded")
            res.ref_uid = None
    return notes


def relink_abundance(results, spectrum_uids, maskset_uids, atom_paths):
    """An abundance result names its reference observation by uid alone."""
    return _relink_reference(results, spectrum_uids, "abundance")


def relink_corrections(results, spectrum_uids, maskset_uids, atom_paths):
    """A correction carries its own reference spectrum and names its atom by path; one
    whose atom is not restored stays listed in the trailing Corrections section."""
    notes = _relink_reference(results, spectrum_uids, "abundcorr")
    for res in results:
        if res.atom_path and res.atom_path not in atom_paths:
            notes.append(f"correction {res.name or res.uid}: model atom "
                         f"{res.atom_path} is not loaded; listed under Corrections")
    return notes


def relink_balance(results, spectrum_uids, maskset_uids, atom_paths):
    """Besides the reference observation, a balance names the mask set it was measured
    from (as a string uid); the link is cleared when that set is not restored, so the panel
    stops offering a refit against windows that are not there."""
    notes = _relink_reference(results, spectrum_uids, "ionbal")
    known = {str(uid) for uid in maskset_uids}
    for res in results:
        if res.maskset_uid and res.maskset_uid not in known:
            notes.append(f"balance {res.name or res.uid}: line mask "
                         f"{res.maskset_name or res.maskset_uid} is not loaded")
            res.maskset_uid = ""
    return notes


# the result-list sections of a saved workspace (workspace.session.RESULT_TOOLS), by the
# name of the tool that owns each
FACTORIES = {"abundance": ab.AbundanceResult.from_dict,
             "abundcorr": ab.AbundanceResult.from_dict,
             "ionbal": ms.BalanceResult.from_dict}
RELINKERS = {"abundance": relink_abundance,
             "abundcorr": relink_corrections,
             "ionbal": relink_balance}


def deserialise(name, records):
    """The result objects of section ``name`` of a saved envelope (no live state touched);
    a malformed record raises."""
    return [FACTORIES[name](d) for d in records]


def relink(name, results, spectrum_uids, maskset_uids, atom_paths):
    """Clear the references of section ``name``'s results that name objects a load does not
    restore; returns one line per cleared reference."""
    return RELINKERS[name](results, spectrum_uids, maskset_uids, atom_paths)
