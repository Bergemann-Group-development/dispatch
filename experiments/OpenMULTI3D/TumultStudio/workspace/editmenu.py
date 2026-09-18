"""The toolbar's Edit menu: **Rename…**, **Duplicate** and **Delete** on the focused object.

The menu acts on the object holding the single left-sidebar focus (``workspace.focus_kind``): a
spectrum, a line mask set, a result or a model atom. Each kind registers a :class:`Target`
saying how its objects are found, renamed and duplicated; an action a kind lacks is shown
disabled. A model atom is identified by its file (loading the file again replaces it), so it
offers Delete only.

The items are rendered here with the focused object's kind, uid and name, and rendered
again whenever the focus moves and when the menu opens. The page asks for a new name with
the browser's prompt and writes ``"<op>\\0<kind>\\0<uid>\\0<name>\\0<nonce>"`` to the
``edit-menu-actions`` channel (NUL-separated, since a name may hold spaces and colons).
Delete takes the path of the sidebar remove buttons (:mod:`removal`), so it asks first when
the removal costs something and offers Undo.
"""
from __future__ import annotations

import html
from dataclasses import dataclass
from typing import Callable

CHANNEL = "edit-menu-actions"
CONTAINER = "edit-object-menu"


@dataclass
class Target:
    """One kind of focusable object. ``find(uid)`` returns the object (with a ``name``) or
    ``None``; ``rename(obj, name)`` and ``duplicate(obj)`` are absent where the kind has no
    such action; ``removal`` is the kind's name in :mod:`removal`."""
    find: Callable
    removal: str
    rename: Callable | None = None
    duplicate: Callable | None = None
    fixed_note: str = ""       # why a disabled Rename / Duplicate is unavailable


ITEMS = (("rename", "Rename…", "rename"), ("duplicate", "Duplicate", "duplicate"),
         ("delete", "Delete", None))


class EditMenu:
    def __init__(self, ctx, focused):
        """``focused()`` returns the ``(kind, uid)`` holding the sidebar focus."""
        self.ctx, self.focused = ctx, focused
        self.targets = {}
        self.bridge = ctx.bridges.html(CONTAINER)
        self.actions = ctx.bridges.action(CHANNEL)
        self.actions.on_change("value", self._on_action)
        ctx.subscribe("focus-changed", lambda **_: self.render())
        ctx.subscribe("masks-changed", lambda **_: self.render())

    def register(self, kind, target):
        self.targets[kind] = target

    def current(self):
        """``(kind, target, obj)`` of the focused object; ``obj`` is ``None`` when nothing
        is focused."""
        kind, uid = self.focused()
        target = self.targets.get(kind)
        obj = target.find(uid) if target is not None and uid is not None else None
        return kind, target, obj

    def html(self):
        kind, target, obj = self.current()
        items = []
        for op, label, attr in ITEMS:
            if obj is None:
                title = "Select an object in the sidebar first"
            elif attr is not None and getattr(target, attr) is None:
                title = target.fixed_note
            else:
                title = ""
            data = ("" if title else
                    f' data-kind="{html.escape(kind)}" data-uid="{html.escape(str(obj.uid))}"'
                    f' data-name="{html.escape(obj.name)}"')
            state = f' disabled title="{html.escape(title)}"' if title else ""
            items.append(f'<button class="menu-item" role="menuitem" data-edit="{op}"'
                         f'{data}{state}>{label}</button>')
        return "".join(items)

    def render(self):
        self.bridge.value = self.html()

    def _on_action(self, attr, old, new):
        if not new:
            return
        op, kind, uid, name = (new.split("\0") + ["", "", "", ""])[:4]
        target = self.targets.get(kind)
        obj = target.find(uid) if target is not None else None
        if obj is not None:
            if op == "rename" and target.rename is not None:
                name = name.strip()
                if name and name != obj.name:
                    target.rename(obj, name)
            elif op == "duplicate" and target.duplicate is not None:
                target.duplicate(obj)
            elif op == "delete":
                self.ctx.removal.request(target.removal, uid)
        self.render()
