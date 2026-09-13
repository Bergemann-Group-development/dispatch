"""The TumultStudio2 config store: truthful read/write outcomes and serialised writers.

The store's temporary directory is the whole world here; nothing touches the project's
own ``.tumultstudio.json``. No physical content is involved: snapshots are plain dicts.
"""
import json
import os
import subprocess
import sys
import threading
from pathlib import Path

import pytest

import store


@pytest.fixture
def config(tmp_path, monkeypatch):
    path = tmp_path / "cfg.json"
    monkeypatch.setattr(store, "CONFIG_PATH", path)
    return path


def _read(path):
    return json.loads(path.read_text())


# --- item 1: missing versus broken config ------------------------------------
def test_missing_config_starts_empty_and_first_write_creates_a_versioned_file(config):
    assert store.list_sessions() == [] and store.recent_files() == []
    assert store.get_pref("x", "dflt") == "dflt"
    slot = store.save_session("first", {"objects": []})
    data = _read(config)
    assert data["version"] == store.CONFIG_VERSION
    assert data["sessions"] == [slot] and slot["id"] == 1
    assert not list(config.parent.glob("*.tmp"))


@pytest.mark.parametrize("text, cause", [
    ("{not json", "not valid JSON"),
    ("[]", "not a JSON object"),
    ('{"preferences": {}, "recent_files": [], "sessions": []}', "version None"),
    ('{"version": 99, "preferences": {}, "recent_files": [], "sessions": []}', "version 99"),
    ('{"version": 2, "preferences": [], "recent_files": [], "sessions": []}',
     "'preferences'"),
    ('{"version": 2, "preferences": {}, "recent_files": [3], "sessions": []}',
     "'recent_files'"),
    ('{"version": 2, "preferences": {}, "recent_files": [], "sessions": [{"id": "1"}]}',
     "malformed slot"),
    ('{"version": 2, "preferences": {}, "recent_files": [], "sessions": ['
     '{"id": 1, "name": "a", "saved_at": 1, "snapshot": {}, "workflows": []},'
     '{"id": 1, "name": "b", "saved_at": 2, "snapshot": {}, "workflows": []}]}', "occurs twice"),
])
def test_broken_config_is_rejected_by_readers_and_writers_alike(config, text, cause):
    """Only an absent file means defaults; anything else names the path and the cause,
    and no writer replaces the broken file with an empty one."""
    config.write_text(text)
    with pytest.raises(store.ConfigError, match=cause) as info:
        store.list_sessions()
    assert str(config) in str(info.value) and info.value.path == config
    for write in (lambda: store.save_session("s", {}), lambda: store.set_pref("k", 1),
                  lambda: store.update_pref("sec", "k", 1), lambda: store.add_recent("f"),
                  lambda: store.delete_session(1)):
        with pytest.raises(store.ConfigError, match=cause):
            write()
    assert config.read_text() == text
    assert not list(config.parent.glob("*.tmp"))


def test_unreadable_config_is_an_error_rather_than_defaults(config):
    config.write_text(json.dumps(store._empty()))
    config.chmod(0)
    try:
        if os.access(config, os.R_OK):
            pytest.skip("cannot make the file unreadable here (running as root?)")
        with pytest.raises(store.ConfigError, match="Could not read config"):
            store.get_pref("x")
        with pytest.raises(store.ConfigError, match="Could not read config"):
            store.set_pref("x", 1)
    finally:
        config.chmod(0o644)
    assert _read(config)["preferences"] == {}


# --- item 2: write and delete outcomes -----------------------------------------
def test_write_failure_names_the_destination_and_keeps_the_previous_file(config, monkeypatch):
    store.save_session("kept", {"n": 1})
    before = config.read_text()

    def fail(src, dst):
        raise OSError(28, "No space left on device")
    monkeypatch.setattr(store.os, "replace", fail)
    with pytest.raises(store.ConfigError, match="No space left on device") as info:
        store.save_session("lost", {"n": 2})
    assert str(config) in str(info.value)
    assert config.read_text() == before
    assert not list(config.parent.glob("*.tmp")), "the writer's temporary file is removed"
    with pytest.raises(store.ConfigError):
        store.delete_session(1)
    assert config.read_text() == before


def test_missing_parent_directory_fails_the_write(tmp_path, monkeypatch):
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "absent" / "cfg.json")
    with pytest.raises(store.ConfigError, match="absent"):
        store.save_session("s", {})
    assert not (tmp_path / "absent").exists()


def test_delete_reports_whether_the_slot_existed(config):
    a = store.save_session("a", {})
    b = store.save_session("b", {})
    assert store.delete_session(a["id"]) is True
    assert [s["id"] for s in store.list_sessions()] == [b["id"]]
    stamp = config.stat().st_mtime_ns
    assert store.delete_session(a["id"]) is False
    assert config.stat().st_mtime_ns == stamp, "a no-op transaction does not rewrite"
    assert store.get_session(a["id"]) is None and store.get_session(b["id"]) == {}


def test_unserialisable_snapshot_is_a_config_error_and_writes_nothing(config):
    store.set_pref("k", 1)
    before = config.read_text()
    with pytest.raises(store.ConfigError, match="not JSON"):
        store.save_session("bad", {"x": object()})
    assert config.read_text() == before


# --- item 3: serialised transactions -----------------------------------------------
def test_threads_saving_slots_and_preferences_lose_nothing(config):
    """Sixteen threads released together each save a slot, touch a shared preferences
    section and add a recent file; every change survives and every slot id is distinct."""
    n = 16
    barrier = threading.Barrier(n, timeout=10)
    errors = []

    def writer(i):
        try:
            barrier.wait()
            store.save_session(f"slot-{i}", {"i": i})
            store.update_pref("shared", f"k{i}", i)
            store.add_recent(f"/f/{i}")
        except Exception as exc:  # noqa: BLE001 - reported through the assertion
            errors.append(exc)

    threads = [threading.Thread(target=writer, args=(i,)) for i in range(n)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(10)
    assert not errors and not any(t.is_alive() for t in threads)
    slots = store.list_sessions()
    assert sorted(s["id"] for s in slots) == list(range(1, n + 1))
    assert {s["name"] for s in slots} == {f"slot-{i}" for i in range(n)}
    assert store.get_section("shared") == {f"k{i}": i for i in range(n)}
    assert len(store.recent_files()) == store.RECENT_MAX
    assert not list(config.parent.glob("*.tmp"))


_PROCESS_WRITER = """
import sys, time
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import store
store.CONFIG_PATH = Path(sys.argv[2])
tag, count, start = sys.argv[3], int(sys.argv[4]), float(sys.argv[5])
time.sleep(max(0.0, start - time.time()))
for i in range(count):
    store.save_session(f"{tag}-{i}", {"tag": tag, "i": i})
    store.update_pref("procs", f"{tag}-{i}", i)
"""


def test_local_processes_sharing_the_config_serialise_their_transactions(config, tmp_path):
    """Four interpreters released at the same instant each save slots and preferences
    against one config path; the file lock keeps every slot and preference."""
    import time
    n_proc, n_each = 4, 10
    script = tmp_path / "writer.py"
    script.write_text(_PROCESS_WRITER)
    start = time.time() + 1.0
    procs = [subprocess.Popen([sys.executable, str(script), str(Path(store.__file__).parent),
                               str(config), f"p{k}", str(n_each), str(start)],
                              stderr=subprocess.PIPE, text=True)
             for k in range(n_proc)]
    for p in procs:
        _, err = p.communicate(timeout=60)
        assert p.returncode == 0, err
    slots = store.list_sessions()
    assert len(slots) == n_proc * n_each
    assert sorted(s["id"] for s in slots) == list(range(1, n_proc * n_each + 1))
    assert len(store.get_section("procs")) == n_proc * n_each
    assert not list(config.parent.glob("*.tmp"))


def test_lock_held_elsewhere_times_out_with_a_clear_error(config, monkeypatch):
    import fcntl
    monkeypatch.setattr(store, "LOCK_TIMEOUT_S", 0.2)
    lock = config.with_name(config.name + ".lock")
    fd = os.open(lock, os.O_RDWR | os.O_CREAT)
    fcntl.flock(fd, fcntl.LOCK_EX)
    try:
        with pytest.raises(store.ConfigError, match="held the config lock"):
            store.set_pref("k", 1)
    finally:
        os.close(fd)
    assert not config.exists()
    store.set_pref("k", 1)  # released: the next transaction goes through
    assert store.get_pref("k") == 1


def test_slots_carry_workflow_references_beside_the_snapshot(config):
    """A slot's workflow references are saved and read back as given, outside the
    snapshot; a malformed reference makes the whole config unreadable."""
    ref = store.workflow_reference("20260901T000000-abc", "/runs/_workflows/x", "t5750")
    assert ref == {"id": "20260901T000000-abc", "folder": "/runs/_workflows/x",
                   "name": "t5750", "imported": False}
    store.save_session("s", {"schema_version": 2}, [ref])
    store.save_session("plain", {"schema_version": 2})
    slot = store.get_slot(1)
    assert slot["workflows"] == [ref] and slot["snapshot"] == {"schema_version": 2}
    assert store.get_slot(2)["workflows"] == [] and store.get_slot(3) is None
    assert store.get_session(1) == {"schema_version": 2}

    data = _read(config)
    data["sessions"][0]["workflows"] = [{"id": 1, "folder": "f", "name": "n", "imported": 0}]
    config.write_text(json.dumps(data))
    with pytest.raises(store.ConfigError, match="malformed workflow reference"):
        store.list_sessions()
    del data["sessions"][0]["workflows"]
    config.write_text(json.dumps(data))
    with pytest.raises(store.ConfigError, match="malformed slot"):
        store.list_sessions()


def test_the_autosave_is_one_slot_with_a_fixed_id_listed_first(config):
    """Autosaves replace one slot with id 0 rather than adding; manual saves keep counting
    from 1 around it; the list puts the autosave ahead of newer manual saves; it loads
    and deletes as any slot."""
    store.save_session("first", {"schema_version": 2, "n": 1})
    auto = store.save_autosave({"schema_version": 2, "n": 2})
    assert auto["id"] == store.AUTOSAVE_ID and auto["autosave"]
    assert auto["name"] == store.AUTOSAVE_NAME
    store.save_autosave({"schema_version": 2, "n": 3}, [store.workflow_reference("w", "/f", "w")])
    store.save_session("second", {"schema_version": 2, "n": 4})
    assert [s["id"] for s in store.list_sessions()] == [0, 2, 1]
    slot = store.autosave_slot()
    assert slot["snapshot"] == {"schema_version": 2, "n": 3} and slot["workflows"][0]["id"] == "w"
    assert store.get_session(0) == slot["snapshot"]
    assert store.delete_session(0) and store.autosave_slot() is None
    assert [s["id"] for s in store.list_sessions()] == [2, 1]
