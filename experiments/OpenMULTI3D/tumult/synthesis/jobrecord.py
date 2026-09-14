"""A dispatch run's record of itself, from two writers with different knowledge.

``dispatch.x`` writes ``job.nml`` into its own run directory as soon as that exists
(``SRC/job_record_mod.f90``): process id, host, working directory, input file, start
time, MPI size, thread count, Slurm job id and code hash, with ``state='running'``; a
normal exit rewrites it as ``finished``. That makes every run discoverable whoever started
it: a shell, a notebook, a batch script or the TumultStudio2 server.

The launcher (``tumult.launch``) writes ``job.json`` beside it with what dispatch cannot
know: that a Slurm job was queued before any process existed (``submitted``), that the
launcher killed the run (``cancelled``), or the exit status of a run that died
(``failed``), plus a digest of the namelist it was launched with.

:func:`read` merges the two into one record and :func:`state_of` resolves its state
against the disk cache, the process table and the scheduler:

* ``finished``   the run wrote ``abund.txt``, or dispatch recorded a normal exit
* ``cancelled`` / ``failed``  the launcher said so
* ``running``    dispatch recorded a start and its process is still alive; also a run
                 recorded on another host, whose process cannot be inspected from here
* ``submitted``  a Slurm job the scheduler still lists (or one it cannot be asked about)
* ``lost``       recorded as running, but the process or Slurm job is gone: a crash,
                 a kill, or a server restarted under it
* ``unknown``    no record and no completed output

Liveness is ``os.kill(pid, 0)`` on the recording host, which cannot distinguish our
process from an unrelated one that reused the id after a reboot. Nothing here acts on the
answer without the further checks the workflow layer makes, so that ambiguity costs a
stale "running" label, not a wrong decision.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import signal
import socket
import subprocess
import time
from pathlib import Path

from tumult.launch import run_complete
from tumult.namelist import NameList

JOB_NAME = "job.json"       # the launcher's record
NATIVE_NAME = "job.nml"     # dispatch's own record


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


def finish(folder, state, **details):
    """Mark a recorded run ended. ``details`` (exit status, log path, error text) are
    stored with the state so a failed run keeps its diagnostics. Absent record, nothing to do."""
    path = Path(folder) / JOB_NAME
    try:
        record = json.loads(path.read_text())
        record["state"] = state
        record["ended"] = time.time()
        record.update(details)
        _write_json(path, record)
    except (OSError, ValueError):
        logging.debug("no %s to finish in %s", JOB_NAME, folder)


def read_launcher(folder):
    try:
        return json.loads((Path(folder) / JOB_NAME).read_text())
    except (OSError, ValueError):
        return None


def read_native(folder):
    """Dispatch's own record, or ``None``. Empty strings become ``None``."""
    path = Path(folder) / NATIVE_NAME
    if not path.is_file():
        return None
    try:
        group = NameList.from_file(path).get("job_record") or {}
    except Exception:       # a record being written, or a namelist this parser rejects
        logging.debug("unreadable %s", path, exc_info=True)
        return None
    native = {key: (None if isinstance(value, str) and not value.strip() else value)
              for key, value in group.items()}
    return native or None


def read(folder):
    """The merged record: the launcher's fields, overlaid with what dispatch recorded
    about itself. ``pid``, ``host`` and ``slurm_id`` come from dispatch when it wrote
    them, since under ``sbatch`` or ``mpirun`` the launcher's child is not the run.
    ``state`` is the raw claim; :func:`state_of` resolves it."""
    launcher = read_launcher(folder)
    native = read_native(folder)
    if launcher is None and native is None:
        return None
    record = dict(launcher or {})
    record["output"] = str(folder)
    if native:
        record["native"] = native
        record["launcher_pid"] = record.get("pid")
        record.update(pid=native.get("pid"), host=native.get("host"),
                      cwd=native.get("cwd"), input=native.get("input"),
                      mpi_size=native.get("mpi_size"), omp_threads=native.get("omp_threads"),
                      dispatch_hash=native.get("dispatch_hash"),
                      slurm_id=native.get("slurm_job_id") or record.get("slurm_id"))
        record.setdefault("started", native.get("started"))
        launched = (launcher or {}).get("state")
        if launched in ("cancelled", "failed"):
            record["state"] = launched      # the launcher saw the end dispatch could not
        else:
            record["state"] = "finished" if native.get("state") == "finished" else "running"
    return record


def process_runs(pid, marker):
    """Whether ``pid`` is alive and its command line still names ``marker``, which
    rules out an unrelated process that reused the id. ``None`` when the command line
    cannot be inspected on this platform."""
    try:
        os.kill(int(pid), 0)
    except (OSError, TypeError, ValueError):
        return False
    if os.name != "posix":
        return None
    try:
        out = subprocess.run(["ps", "-o", "command=", "-p", str(int(pid))],
                             capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    return marker in out


def cancel(folder, record=None):
    """End a recorded run that is still going: ``scancel`` its Slurm job, else SIGTERM
    its process on this host once the process table confirms the pid still runs with the
    recorded input file on its command line. Returns whether anything was signalled; the
    launcher record is then marked ``cancelled`` so the folder does not read as lost."""
    record = record if record is not None else read(folder)
    if record is None or state_of(folder, record) not in ("running", "submitted"):
        return False
    if record.get("slurm_id"):
        from . import scheduler
        try:
            sent = scheduler.scancel(record["slurm_id"])
        except OSError:
            return False
    else:
        pid, marker = record.get("pid"), Path(record.get("input") or "").name
        if not pid or not marker or record.get("host") not in (None, socket.gethostname()):
            return False
        if process_runs(pid, marker) is not True:
            return False
        try:
            os.kill(int(pid), signal.SIGTERM)
        except OSError:
            return False
        sent = True
    if sent:
        if read_launcher(folder) is None:
            write(folder, {}, pid=record.get("pid"), slurm_id=record.get("slurm_id"))
        finish(folder, "cancelled")
    return sent


def alive(pid, host=None):
    """Whether ``pid`` is alive on this host; ``None`` when the record names another
    host, where the process table cannot be consulted from here."""
    if not pid:
        return False
    if host and host != socket.gethostname():
        return None
    try:
        os.kill(int(pid), 0)
    except (OSError, TypeError, ValueError):
        return False
    return True


def state_of(folder, record=None):
    """Resolve a record's state against the disk cache, the process table and the
    scheduler (see the module docstring)."""
    record = record if record is not None else read(folder)
    if record is None:
        return "finished" if run_complete(folder) else "unknown"
    if run_complete(folder):
        return "finished"
    state = record.get("state")
    if state not in ("running", "submitted"):
        return state
    if record.get("slurm_id"):
        from . import scheduler
        try:
            listed = scheduler.state(record["slurm_id"])
        except OSError:
            return "submitted"      # no scheduler client here; nothing to poll
        if listed is None:
            return "lost"
        return "running" if listed == "running" and state == "running" else "submitted"
    if state == "submitted":
        return "submitted"
    live = alive(record.get("pid"), record.get("host"))
    if live is None:
        return "running"            # another host: trust the record
    return "running" if live else "lost"


def scan(root, limit=None):
    """Every recorded job under ``root``, newest first, with its resolved state.
    Runs with only dispatch's own record (started outside the launcher) are included."""
    found = {}
    root = Path(root)
    for name in (JOB_NAME, NATIVE_NAME):
        for path in root.rglob(name):
            if path.parent in found:
                continue
            record = read(path.parent)
            if record is None:
                continue
            record["resolved"] = state_of(path.parent, record)
            found[path.parent] = record
    records = list(found.values())
    records.sort(key=started_time, reverse=True)
    return records[:limit] if limit else records


def run_id(folder):
    """An opaque, colon-free id for the run in ``folder``, stable across processes, for
    the daemon's job list and the GUI's action channel."""
    return "run-" + hashlib.sha1(str(folder).encode()).hexdigest()[:12]


def started_time(record):
    """A sortable start time (epoch seconds) from either writer: the launcher's epoch
    seconds or dispatch's ISO timestamp; 0 when neither recorded one."""
    value = record.get("started")
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        from datetime import datetime
        try:
            return datetime.fromisoformat(value).timestamp()
        except ValueError:
            return 0.0
    return 0.0
