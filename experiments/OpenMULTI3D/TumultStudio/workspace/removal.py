"""Removing objects from the workspace: the one path every sidebar remove button takes.

A remove button sends ``(kind, uid)`` to :meth:`Removal.request`, which works out what the
removal would cost before doing it:

* *dependents*: results that name the object (a spectrum as their target, a line mask set
  as a balance's windows) and cannot show or recompute their measurement without it, or, for a
  model atom, the correction results that are removed with it;
* *unsaved edits*: a line mask set edited since its file was read or written, or one that
  was never written;
* *running work*: the workflow computing a result, which the removal cancels.

A removal that costs none of these happens at once. Otherwise the page asks first, listing
what is affected (``htmlbridge.confirm_bridge``); the answer comes back as
``confirm <token>`` and is acted on only while the plan is still the one the user saw.

Every removal can be undone from its status message while it is the latest one and the
workspace has not been replaced since: the object returns with its uid, list position,
selection and state, so results that name it work again. That is what keeps the autosave,
which records the removal within seconds, from turning a slip into a loss. The undo does not
restart a cancelled workflow, which is why cancelling one is always confirmed. Only the
latest removal is kept; this is not an undo history.
"""
from __future__ import annotations

import itertools
import json
from dataclasses import dataclass, field
from typing import Callable

from ..ui.htmlbridge import parse_action

CHANNEL = "removal-actions"
UNDO_MS = 10000          # how long the status message offering Undo stays


@dataclass
class Kind:
    """How one kind of object is found, assessed and taken out.

    ``detach(obj)`` removes the object and returns a callable that puts it back."""
    label: str
    find: Callable
    detach: Callable
    dependents: Callable = lambda obj: []
    dependents_note: str = ("These results refer to it and cannot show or recompute "
                            "their measurement while it is gone:")
    unsaved: Callable = lambda obj: ""
    running: Callable = lambda obj: []


@dataclass
class Plan:
    """What removing one object would cost."""
    kind: str
    uid: str
    name: str
    dependents: list = field(default_factory=list)
    dependents_note: str = ""
    unsaved: str = ""
    running: list = field(default_factory=list)

    @property
    def needs_confirmation(self):
        return bool(self.dependents or self.unsaved or self.running)

    def message(self, label):
        parts = [f"Remove the {label} “{self.name}”?"]
        if self.unsaved:
            parts.append(self.unsaved)
        if self.dependents:
            parts.append(self.dependents_note + "\n"
                         + "\n".join(f"  • {n}" for n in self.dependents))
        if self.running:
            parts.append("This calculation is cancelled and cannot be resumed by Undo:\n"
                         + "\n".join(f"  • {n}" for n in self.running))
        parts.append("Undo in the status bar puts it back.")
        return "\n\n".join(parts)


class Removal:
    def __init__(self, ctx, notify):
        self.ctx, self.notify = ctx, notify
        self.kinds = {}
        self._tokens = itertools.count(1)
        self.pending = {}          # token -> (kind, uid, message) awaiting the user's answer
        self._undo = None          # (token, generation, name, restore) of the latest removal
        self.confirm_bridge = ctx.bridges.confirm(CHANNEL)
        self.actions = ctx.bridges.action(CHANNEL)
        self.actions.on_change("value", self._on_action)

    def register(self, name, kind):
        self.kinds[name] = kind

    def _on_action(self, attr, old, new):
        if not new:
            return
        action, arg = parse_action(new)
        if action == "confirm":
            self.confirm(arg)
        elif action == "undo":
            self.undo(arg)

    # --- planning -------------------------------------------------------------------------
    def plan(self, kind, uid):
        """The :class:`Plan` for removing ``uid``, or ``None`` when it is not listed."""
        spec = self.kinds[kind]
        obj = spec.find(uid)
        if obj is None:
            return None
        return Plan(kind, uid, obj.name, dependents=list(spec.dependents(obj)),
                    dependents_note=spec.dependents_note, unsaved=spec.unsaved(obj),
                    running=list(spec.running(obj)))

    def request(self, kind, uid):
        """A remove button: remove at once, or ask first when the removal costs something."""
        plan = self.plan(kind, uid)
        if plan is None:
            return
        if not plan.needs_confirmation:
            self.remove(kind, uid)
            return
        token = str(next(self._tokens))
        message = plan.message(self.kinds[kind].label)
        self.pending = {token: (kind, uid, message)}     # a newer question replaces an older
        self.confirm_bridge.value = (json.dumps([message, f"confirm {token}"])
                                     + f"\0{token}")

    def confirm(self, token):
        """The user accepted question ``token``. The removal goes ahead if what it costs is
        still what the question listed; otherwise the new plan is asked about."""
        asked = self.pending.pop(str(token), None)
        if asked is None:
            return
        kind, uid, message = asked
        plan = self.plan(kind, uid)
        if plan is None:
            return
        if plan.needs_confirmation and plan.message(self.kinds[kind].label) != message:
            self.request(kind, uid)
            return
        self.remove(kind, uid)

    # --- removal and its undo -------------------------------------------------------------
    def remove(self, kind, uid):
        """Remove ``uid`` without asking and offer to undo it."""
        spec = self.kinds[kind]
        obj = spec.find(uid)
        if obj is None:
            return
        name = obj.name
        restore = spec.detach(obj)
        token = str(next(self._tokens))

        def put_back():
            if spec.find(uid) is not None:
                raise ValueError(f"{uid} is listed again")
            restore()
        self._undo = (token, self.ctx.generation, name, put_back)
        self.notify(f"Removed “{name}”.", level="info", revert_ms=UNDO_MS,
                    undo=token)

    def undo(self, token):
        """Put back the removal ``token`` if it is the latest and the workspace it came from
        is still the live one."""
        if self._undo is None or self._undo[0] != str(token):
            return
        _, generation, name, put_back = self._undo
        self._undo = None
        if generation is not self.ctx.generation:
            self.notify(f"“{name}” belonged to a workspace that has since been "
                        "replaced; it cannot be put back.")
            return
        try:
            put_back()
        except ValueError:
            self.notify(f"“{name}” cannot be put back: an object with its "
                        "identity is listed again.")
            return
        self.notify(f"Restored “{name}”.", level="success")


def register_workspace(removal, canvas, masktool, wizards, atomtool):
    """The workspace's removable kinds: spectra, line mask sets, the analysis results
    (``wizards`` by tool name) and model atoms."""
    def results_naming(test):
        return [r.name for w in wizards.values() for r in w.results if test(r)]

    def detach_spectrum(obj):
        index = canvas.objects.index(obj)
        selected = canvas.selected_uid == obj.uid
        canvas.remove_object(obj.uid)
        return lambda: canvas.restore_object(obj, index, select=selected)

    removal.register("spectrum", Kind(
        label="spectrum", find=canvas._find, detach=detach_spectrum,
        dependents=lambda obj: results_naming(lambda r: r.ref_uid == obj.uid)))

    def find_set(uid):
        try:
            return masktool.store.find(int(uid))
        except (TypeError, ValueError):
            return None

    def unsaved_set(mset):
        if mset.dirty:
            return "Its edits have not been saved to its file."
        if not mset.baseline and mset.masks:
            return "It has never been saved to a file."
        return ""

    def detach_set(mset):
        placed = masktool.remove_set(mset.uid)
        return lambda: masktool.store.insert_set(**placed)

    removal.register("linemask", Kind(
        label="line mask set", find=find_set, detach=detach_set, unsaved=unsaved_set,
        # only a balance reads its windows back from the set (IonBalanceTool._source_set)
        dependents=lambda mset: [r.name for r in wizards["ionbal"].results
                                 if r.maskset_uid == str(mset.uid)]))

    for name, wizard in wizards.items():
        def detach_result(res, wizard=wizard):
            placed = wizard.remove_result(res.uid)
            return lambda: wizard.insert_result(**placed)

        removal.register(name, Kind(
            label="result", find=wizard._find, detach=detach_result,
            running=lambda res, wizard=wizard: [res.name] if wizard.is_computing(res) else []))

    corr = atomtool.corr

    def detach_atom(atom):
        placed = atomtool.remove(atom.uid)
        return lambda: atomtool.insert(**placed)

    removal.register("atom", Kind(
        label="model atom", find=atomtool.find, detach=detach_atom,
        dependents=lambda atom: [r.name for r in corr.results_for(atom.path)],
        dependents_note="These correction results were computed from it and are removed "
                        "with it:",
        running=lambda atom: [r.name for r in corr.results_for(atom.path)
                              if corr.is_computing(r)]))
