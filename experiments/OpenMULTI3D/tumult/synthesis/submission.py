"""Durable asynchronous submission of the same resolved calculation used synchronously."""
from __future__ import annotations

from concurrent.futures import Future
from copy import deepcopy
import json
import os
from pathlib import Path
import threading
import time
import uuid

from .cancel import CancelToken, CancelledError
from .request import ResolvedSynthRequest, SynthesisResult, SynthRequest
from .runtime import WORKFLOWS, observe_runs

RECORD = 'workflow.json'
_live = {}
_lock = threading.RLock()


def _atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'{path.name}.{os.getpid()}.partial')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + '\n')
    tmp.replace(path)


def _read(folder):
    try:
        return json.loads((Path(folder) / RECORD).read_text())
    except (OSError, ValueError):
        return None


class WorkflowHandle:
    """Stable job identity with status, cancellation and durable result retrieval."""

    def __init__(self, folder, future=None, token=None):
        self.folder = Path(folder).expanduser().resolve()
        self.id = self.folder.name
        self._future = future
        self._token = token

    def _record(self):
        record = _read(self.folder)
        if record is None:
            raise FileNotFoundError(f'No workflow record in {self.folder}')
        if (self.folder / 'result' / 'result.json').is_file():
            record['state'] = 'finished'
        elif record.get('state') in ('queued', 'running') and self._future is None:
            record['state'] = 'interrupted'
        return record

    @property
    def status(self):
        return self._record()['state']

    def _native_runs(self, record):
        """State and progress of each native run this workflow has scheduled so far."""
        from . import jobrecord
        from tumult.launch import monitor
        items = []
        for folder in record.get('runs') or ():
            folder = Path(folder)
            log = folder.with_suffix('.log')
            native = jobrecord.read(folder)
            items.append({'output': str(folder), 'log': str(log) if log.is_file() else None,
                          'state': jobrecord.state_of(folder, native),
                          'progress': monitor(folder, log if log.is_file() else None)})
        return items

    @property
    def progress(self):
        record = self._record()
        diagnostics = deepcopy(record.get('diagnostics') or {})
        phase = record.get('phase')
        if record.get('state') == 'running':
            try:
                native = self._native_runs(record)
            except (OSError, KeyError, ValueError):
                native = []
            if native:
                diagnostics.update(
                    native_runs=len(native),
                    completed=sum(item['state'] == 'finished' for item in native),
                    done=sum(item['progress'][0] for item in native),
                    total=sum(item['progress'][1] for item in native),
                    estimated=True, runs=native)
        return {'phase': phase, 'elapsed': max(
            0.0, (record.get('ended') or time.time()) - record['submitted']),
                'diagnostics': diagnostics}

    @property
    def runs(self):
        return tuple(self._record().get('runs') or ())

    @property
    def logs(self):
        logs = [str(Path(run).with_suffix('.log')) for run in self.runs
                if Path(run).with_suffix('.log').is_file()]
        return tuple(dict.fromkeys(logs))

    def cancel(self):
        """Cancel this subscriber's calculation; never signal an unverified saved PID."""
        if self._token is None or self._future is None or self._future.done():
            return False
        self._token.cancel()
        return True

    def result(self, timeout=None):
        if self._future is not None:
            self._future.result(timeout=timeout)
        state = self.status
        if state != 'finished':
            raise RuntimeError(f'Workflow {self.id} is {state}')
        return SynthesisResult.load(self.folder / 'result')

    def to_dict(self):
        return self._record()


def submit(request, *, root=None):
    """Copy, resolve and queue ``request``; return immediately with a durable handle."""
    from .workflows import execute_request, resolve_request
    resolved = request if isinstance(request, ResolvedSynthRequest) else resolve_request(request)
    # JSON round-trip proves that no mutable masks/settings remain shared with the caller.
    owned = SynthRequest.from_dict(resolved.request.to_dict())
    resolved = resolve_request(owned)
    base = Path(root or resolved.output_root) / '_workflows'
    identity = f'{time.strftime("%Y%m%dT%H%M%S")}-{uuid.uuid4().hex[:10]}'
    folder = base / identity
    record = {'version': 1, 'id': identity, 'state': 'queued', 'phase': 'queued',
              'submitted': time.time(), 'request': resolved.to_dict(), 'runs': [],
              'result': None}
    _atomic(folder / RECORD, record)
    token = CancelToken()
    handle = WorkflowHandle(folder, token=token)

    def update(**values):
        with _lock:
            current = _read(folder) or record
            current.update(values)
            _atomic(folder / RECORD, current)

    runs = []

    def note(output):
        # Called by the batch executor before each native run starts, cache hits included,
        # so the record names its runs while they are in flight.
        if output not in runs:
            runs.append(output)
            update(runs=list(runs))

    def work():
        update(state='running', phase='preparation', started=time.time())
        try:
            token.check()
            update(phase='synthesis')
            with observe_runs(note):
                result = execute_request(resolved, cancel_token=token)
            update(phase='saving')
            result_path = result.save(folder / 'result')
            update(state='finished', phase='finished', ended=time.time(),
                   result=str(result_path), runs=list(dict.fromkeys(runs + list(result.runs))))
            return str(result_path)
        except CancelledError:
            update(state='cancelled', phase='cancelled', ended=time.time())
            raise
        except BaseException as exc:
            update(state='failed', phase='failed', ended=time.time(),
                   error=f'{type(exc).__name__}: {exc}')
            raise
        finally:
            with _lock:
                _live.pop(identity, None)

    with _lock:
        _live[identity] = handle
        handle._future = WORKFLOWS.submit(work)
    return handle


def open_workflow(path, *, root=None):
    """Open a workflow by folder/path or stable id after a UI reconnect."""
    path = Path(path)
    folder = ((Path(root) / '_workflows' / path) if root is not None and len(path.parts) == 1
              else path).expanduser().resolve()
    with _lock:
        return _live.get(folder.name) or WorkflowHandle(folder)


def discover_workflows(root, limit=None):
    """Newest workflow handles, including recoverable completed results."""
    folders = [path.parent for path in Path(root).expanduser().resolve().rglob(RECORD)]
    folders.sort(key=lambda path: (_read(path) or {}).get('submitted', 0), reverse=True)
    handles = [open_workflow(path) for path in folders]
    return handles[:limit] if limit else handles


def run(request, *, cancel_token=None):
    """Blocking counterpart to :func:`submit`."""
    from .workflows import execute_request
    return execute_request(request, cancel_token=cancel_token)
