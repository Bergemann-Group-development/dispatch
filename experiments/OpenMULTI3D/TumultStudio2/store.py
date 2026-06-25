"""Project-local persistence.

App preferences, recent files, and saved *sessions* all live in
``./.tumultstudio.json`` (override with the ``TUMULT_CONFIG`` env var); writes are
atomic and every accessor reads fresh, so concurrent Bokeh sessions stay consistent.

Sessions are named save *slots* (think computer-game saves) kept inside the config
file itself — never as individual files. Each slot is ``{"id", "name", "saved_at",
"snapshot"}``; manage them via list_sessions/save_session/delete_session.
"""
from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path

CONFIG_PATH = Path(os.environ.get("TUMULT_CONFIG", Path.cwd() / ".tumultstudio.json"))
RECENT_MAX = 15
_SECTIONS = {"preferences": dict, "recent_files": list, "sessions": list}


def _load() -> dict:
    try:
        data = json.loads(CONFIG_PATH.read_text())
    except (OSError, ValueError):
        data = {}
    for key, kind in _SECTIONS.items():
        if not isinstance(data.get(key), kind):
            data[key] = kind()
    return data


def _save(data: dict) -> None:
    tmp = CONFIG_PATH.with_name(CONFIG_PATH.name + ".tmp")
    try:
        tmp.write_text(json.dumps(data, indent=2))
        tmp.replace(CONFIG_PATH)
    except OSError:
        logging.exception("Could not write config %s", CONFIG_PATH)


# --- preferences ------------------------------------------------------------
def get_pref(key, default=None):
    return _load()["preferences"].get(key, default)


def set_pref(key, value) -> None:
    data = _load()
    data["preferences"][key] = value
    _save(data)


def get_section(name) -> dict:
    """A preferences sub-dict as a fresh dict (never ``None``) — handy to merge over defaults."""
    val = _load()["preferences"].get(name)
    return dict(val) if isinstance(val, dict) else {}


def update_pref(section, key, value) -> None:
    """Set one ``key`` inside a preferences sub-dict (atomic read-modify-write)."""
    data = _load()
    sec = data["preferences"].get(section)
    if not isinstance(sec, dict):
        sec = {}
    sec[key] = value
    data["preferences"][section] = sec
    _save(data)


def set_section(name, values: dict) -> None:
    """Replace a whole preferences sub-dict (for components that persist their full state)."""
    data = _load()
    data["preferences"][name] = dict(values)
    _save(data)


# --- recent files (most-recent-first) ---------------------------------------
def recent_files() -> list[str]:
    return _load()["recent_files"]


def add_recent(path) -> None:
    p = str(Path(path))
    data = _load()
    recents = [r for r in data["recent_files"] if r != p]
    recents.insert(0, p)
    data["recent_files"] = recents[:RECENT_MAX]
    _save(data)


# --- sessions (named save slots, stored in the config file) -----------------
def list_sessions() -> list[dict]:
    """Saved slots, most-recently-saved first."""
    slots = _load()["sessions"]
    return sorted(slots, key=lambda s: s.get("saved_at", 0), reverse=True)


def save_session(name, snapshot) -> dict:
    """Add a new save slot and return it. Names need not be unique."""
    data = _load()
    slots = data["sessions"]
    sid = max((s.get("id", 0) for s in slots), default=0) + 1
    slot = {"id": sid, "name": str(name), "saved_at": time.time(), "snapshot": snapshot}
    slots.append(slot)
    data["sessions"] = slots
    _save(data)
    return slot


def get_session(sid) -> dict | None:
    """The snapshot for slot ``sid`` (or ``None`` if it's gone)."""
    for s in _load()["sessions"]:
        if s.get("id") == sid:
            return s.get("snapshot")
    return None


def delete_session(sid) -> None:
    data = _load()
    data["sessions"] = [s for s in data["sessions"] if s.get("id") != sid]
    _save(data)
