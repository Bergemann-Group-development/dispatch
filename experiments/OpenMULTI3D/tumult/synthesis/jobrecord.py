"""A dispatch run's record of itself, written when it starts rather than when it finishes.

The synthesis cache is durable but result-only: a run appears on disk once ``abund.txt``
exists, and until then nothing records that it was ever started. That is adequate for the
sub-minute jobs the GUI runs today, and not for a 3D NLTE departure solve, where a user who
reopens the browser needs to be told "still running" rather than "nothing here".

This writes ``job.json`` into the run's own output folder at launch time and updates it
when the run ends. ``scan`` walks a cache root and reports what it finds, resolving each
record's state against the disk cache and the operating system rather than trusting what
the record last claimed:

* ``finished``  the run wrote ``abund.txt``; it is a cache hit whatever the record says
* ``running``   the record says running and its process is still alive
* ``lost``      the record says running and its process is gone: killed, or the server
                was restarted under it
* ``failed`` / ``cancelled``  the launcher said so before it exited

Liveness is ``os.kill(pid, 0)``, which cannot distinguish our process from an unrelated one
that reused the id after a reboot. Nothing here acts on the answer — it is displayed — so
that ambiguity costs a stale "running" label, not a wrong decision.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from pathlib import Path

from tumult.launch import run_complete

JOB_NAME = "job.json"


def _write_json(path, value):
    path = Path(path)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.partial")
    tmp.write_text(json.dumps(value, indent=1, sort_keys=True) + "\n")
    tmp.replace(path)


def settings_hash(namelist):
    """A stable digest of the namelist a run was launched with."""
    return hashlib.sha256(
        json.dumps(namelist, sort_keys=True, default=str).encode()).hexdigest()[:16]


def write(folder, namelist, *, pid=None, slurm_id=None, state="running"):
    """Record a starting run in its own output folder. Never raises."""
    folder = Path(folder)
    record = {"state": state, "started": time.time(), "pid": pid, "slurm_id": slurm_id,
              "output": str(folder), "settings": settings_hash(namelist)}
    try:
        folder.mkdir(parents=True, exist_ok=True)
        _write_json(folder / JOB_NAME, record)
    except OSError:
        # The record is an observability convenience; a read-only or vanished output
        # directory must not fail the run it describes.
        logging.exception("could not write %s in %s", JOB_NAME, folder)
    return record


def finish(folder, state):
    """Mark a recorded run ended. Absent record, nothing to do."""
    path = Path(folder) / JOB_NAME
    try:
        record = json.loads(path.read_text())
        record["state"] = state
        record["ended"] = time.time()
        _write_json(path, record)
    except (OSError, ValueError):
        logging.debug("no %s to finish in %s", JOB_NAME, folder)


def read(folder):
    try:
        return json.loads((Path(folder) / JOB_NAME).read_text())
    except (OSError, ValueError):
        return None


def alive(pid):
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except (OSError, TypeError, ValueError):
        return False
    return True


def state_of(folder, record=None):
    """Resolve a record's state against the disk cache and the process table."""
    record = record if record is not None else read(folder)
    if record is None:
        return "finished" if run_complete(folder) else "unknown"
    if run_complete(folder):
        return "finished"
    if record.get("state") != "running":
        return record["state"]
    if record.get("slurm_id"):
        return "submitted"      # the scheduler owns it; this process cannot poll it
    return "running" if alive(record.get("pid")) else "lost"


def scan(root, limit=None):
    """Every recorded job under ``root``, newest first, with its resolved state."""
    found = []
    for path in Path(root).rglob(JOB_NAME):
        record = read(path.parent)
        if record is None:
            continue
        record["output"] = str(path.parent)
        record["resolved"] = state_of(path.parent, record)
        found.append(record)
    found.sort(key=lambda r: r.get("started", 0), reverse=True)
    return found[:limit] if limit else found
