"""Project-local persistence.

App preferences, recent files and saved workspace *slots* all live in one JSON config
file, ``./.tumultstudio.json`` by default (override with the ``TUMULT_CONFIG`` env var).

Every writer runs one serialised read-modify-write transaction (`transaction`): a
process-wide lock plus an advisory ``flock`` on ``<config>.lock`` covers threads and
other local server processes sharing the path, so slot IDs are allocated under the lock
and no writer overwrites another's change. The file is written to a writer-specific
temporary file and atomically replaced; a transaction that cannot read the existing
config never writes. Readers take no lock: replacement is atomic, so a read sees a
complete old or new file.

Config contract: a missing file means empty defaults. An unreadable, malformed or
structurally invalid file, or one with an unsupported ``version``, raises `ConfigError`
naming the path and the cause; nothing starts from defaults in that case. Slots are
``{"id", "name", "saved_at", "snapshot", "workflows"}`` inside the file; the snapshot
content is the workspace schema's business (`session.py`) and is not validated here.
``workflows`` lists the durable synthesis workflows the saved document was attached to,
beside the workspace rather than inside it: ``{"id", "folder", "name", "imported"}``
each, where ``imported`` records that the workflow's result already reached the
workspace, so a load offers it for import only once.
"""
from __future__ import annotations

import contextlib
import fcntl
import json
import os
import threading
import time
from pathlib import Path

CONFIG_PATH = Path(os.environ.get("TUMULT_CONFIG", Path.cwd() / ".tumultstudio.json"))
CONFIG_VERSION = 2
RECENT_MAX = 15
LOCK_TIMEOUT_S = 10.0
_SECTIONS = {"preferences": dict, "recent_files": list, "sessions": list}
_SLOT_FIELDS = {"id": int, "name": str, "saved_at": (int, float), "workflows": list}
AUTOSAVE_NAME = "Last session"      # the one slot the running document keeps current
AUTOSAVE_ID = 0                     # its fixed id; manual saves count up from 1
_WORKFLOW_FIELDS = {"id": str, "folder": str, "name": str, "imported": bool}
_thread_lock = threading.Lock()


class ConfigError(Exception):
    """A config read, write, delete or lock failed; the message names the path and cause."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"{message} ({self.path})")


def _empty() -> dict:
    return {"version": CONFIG_VERSION, **{key: kind() for key, kind in _SECTIONS.items()}}


def _validate(data, path) -> dict:
    """The parsed config if it has the supported shape; `ConfigError` otherwise."""
    if not isinstance(data, dict):
        raise ConfigError(path, "Config is not a JSON object")
    version = data.get("version")
    if version != CONFIG_VERSION:
        raise ConfigError(path, f"Config version {version!r} is not supported "
                                f"(expected {CONFIG_VERSION}); move the file aside "
                                f"to start afresh")
    for key, kind in _SECTIONS.items():
        if not isinstance(data.get(key), kind):
            raise ConfigError(path, f"Config section {key!r} is missing or not a "
                                    f"{kind.__name__}")
    if not all(isinstance(r, str) for r in data["recent_files"]):
        raise ConfigError(path, "Config section 'recent_files' holds a non-string entry")
    ids = set()
    for slot in data["sessions"]:
        if not isinstance(slot, dict) or "snapshot" not in slot or any(
                not isinstance(slot.get(field), kind) or isinstance(slot.get(field), bool)
                for field, kind in _SLOT_FIELDS.items()):
            raise ConfigError(path, "Config section 'sessions' holds a malformed slot")
        if not isinstance(slot.get("autosave", False), bool):
            raise ConfigError(path, f"Config slot {slot['id']} has a malformed autosave flag")
        if slot["id"] in ids:
            raise ConfigError(path, f"Config slot id {slot['id']} occurs twice")
        ids.add(slot["id"])
        for ref in slot["workflows"]:
            if not isinstance(ref, dict) or any(
                    not isinstance(ref.get(field), kind)
                    or (kind is not bool and isinstance(ref.get(field), bool))
                    for field, kind in _WORKFLOW_FIELDS.items()):
                raise ConfigError(path, f"Config slot {slot['id']} holds a malformed "
                                        f"workflow reference")
    return data


def _load() -> dict:
    """The config as a validated dict; empty defaults only when the file is absent."""
    path = CONFIG_PATH
    try:
        text = path.read_text()
    except FileNotFoundError:
        return _empty()
    except OSError as exc:
        raise ConfigError(path, f"Could not read config: {exc.strerror or exc}") from exc
    try:
        data = json.loads(text)
    except ValueError as exc:
        raise ConfigError(path, f"Config is not valid JSON: {exc}") from exc
    return _validate(data, path)


def _write(text: str) -> None:
    """Replace the config atomically through a writer-specific temporary file."""
    path = CONFIG_PATH
    tmp = path.with_name(f"{path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    try:
        with open(tmp, "w") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except OSError as exc:
        with contextlib.suppress(OSError):
            tmp.unlink()
        raise ConfigError(path, f"Could not write config: {exc.strerror or exc}") from exc


@contextlib.contextmanager
def _file_lock():
    lock_path = CONFIG_PATH.with_name(CONFIG_PATH.name + ".lock")
    try:
        fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o644)
    except OSError as exc:
        raise ConfigError(CONFIG_PATH, f"Could not open the config lock {lock_path}: "
                                       f"{exc.strerror or exc}") from exc
    try:
        deadline = time.monotonic() + LOCK_TIMEOUT_S
        while True:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise ConfigError(CONFIG_PATH, f"Another writer has held the config "
                                                   f"lock for over {LOCK_TIMEOUT_S:g} s")
                time.sleep(0.01)
        yield
    finally:
        os.close(fd)  # releases the flock


def _dump(data: dict) -> str:
    return json.dumps(data, indent=2)


@contextlib.contextmanager
def transaction():
    """Yield the validated config for one serialised read-modify-write and replace the
    file on normal exit if the dict changed. Raises `ConfigError` when the config cannot
    be read, locked or written; a read failure means nothing is written. Not reentrant."""
    with _thread_lock, _file_lock():
        data = _load()
        before = _dump(data)
        yield data
        try:
            after = _dump(data)
        except (TypeError, ValueError) as exc:
            raise ConfigError(CONFIG_PATH, f"Config content is not JSON: {exc}") from exc
        if after != before:
            _write(after)


# --- preferences ------------------------------------------------------------
def get_pref(key, default=None):
    return _load()["preferences"].get(key, default)


def set_pref(key, value) -> None:
    with transaction() as data:
        data["preferences"][key] = value


def get_section(name) -> dict:
    """A preferences sub-dict as a fresh dict (never ``None``), handy to merge over defaults."""
    val = _load()["preferences"].get(name)
    return dict(val) if isinstance(val, dict) else {}


def update_pref(section, key, value) -> None:
    """Set one ``key`` inside a preferences sub-dict."""
    with transaction() as data:
        sec = data["preferences"].get(section)
        if not isinstance(sec, dict):
            sec = {}
        sec[key] = value
        data["preferences"][section] = sec


def set_section(name, values: dict) -> None:
    """Replace a whole preferences sub-dict (for components that persist their full state)."""
    with transaction() as data:
        data["preferences"][name] = dict(values)


# --- recent files (most-recent-first): the spectrum picker's history ---------
def recent_files() -> list[str]:
    return _load()["recent_files"]


def add_recent(path) -> None:
    p = str(Path(path))
    with transaction() as data:
        recents = [r for r in data["recent_files"] if r != p]
        recents.insert(0, p)
        data["recent_files"] = recents[:RECENT_MAX]


# --- sessions (named save slots, stored in the config file) -----------------
def list_sessions() -> list[dict]:
    """Saved slots, most-recently-saved first, the autosave slot ahead of all of them."""
    return sorted(_load()["sessions"],
                  key=lambda s: (s.get("autosave", False), s["saved_at"]), reverse=True)


def workflow_reference(id, folder, name, imported=False) -> dict:
    """One slot entry naming a durable workflow the document is attached to."""
    return {"id": str(id), "folder": str(folder), "name": str(name), "imported": bool(imported)}


def save_session(name, snapshot, workflows=()) -> dict:
    """Add a new save slot and return it once the config has been replaced on disk.
    Names need not be unique; the id is allocated inside the transaction. ``workflows``
    are `workflow_reference` entries saved beside the snapshot."""
    with transaction() as data:
        slots = data["sessions"]
        sid = max((s["id"] for s in slots), default=AUTOSAVE_ID) + 1
        slot = {"id": sid, "name": str(name), "saved_at": time.time(), "snapshot": snapshot,
                "workflows": [dict(ref) for ref in workflows]}
        slots.append(slot)
    return slot


def save_autosave(snapshot, workflows=()) -> dict:
    """Replace the autosave slot with ``snapshot`` (creating it the first time) and return
    it. It is one slot under :data:`AUTOSAVE_NAME` with the fixed id :data:`AUTOSAVE_ID`,
    flagged ``autosave``, outside the ids manual saves allocate; loading and deleting it
    work as for any other slot."""
    with transaction() as data:
        slots = data["sessions"]
        slot = next((s for s in slots if s.get("autosave")), None)
        if slot is None:
            slot = {"id": AUTOSAVE_ID, "name": AUTOSAVE_NAME, "autosave": True}
            slots.append(slot)
        slot.update(saved_at=time.time(), snapshot=snapshot,
                    workflows=[dict(ref) for ref in workflows])
    return slot


def autosave_slot() -> dict | None:
    """The autosave slot, or ``None`` before the first autosave (or after its deletion)."""
    return next((s for s in _load()["sessions"] if s.get("autosave")), None)


def get_slot(sid) -> dict | None:
    """The whole slot ``sid`` (name, snapshot and workflow references), or ``None``."""
    for s in _load()["sessions"]:
        if s["id"] == sid:
            return s
    return None


def get_session(sid) -> dict | None:
    """The snapshot for slot ``sid`` (or ``None`` if it's gone)."""
    slot = get_slot(sid)
    return None if slot is None else slot["snapshot"]


def delete_session(sid) -> bool:
    """Remove slot ``sid``; whether it existed. The file is untouched when it did not."""
    with transaction() as data:
        before = len(data["sessions"])
        data["sessions"] = [s for s in data["sessions"] if s["id"] != sid]
        return len(data["sessions"]) < before
