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
from .runtime import WORKFLOWS

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

    @property
    def progress(self):
        record = self._record()
        diagnostics = deepcopy(record.get('diagnostics') or {})
        phase = record.get('phase')
        if record.get('state') == 'running':
            try:
                from . import jobrecord
                from tumult.launch import monitor
                root = record['request']['output_root']
                native = [item for item in jobrecord.scan(root)
                          if item.get('started', 0) >= record['submitted'] - 1]
                counts = []
                for item in native:
                    output = Path(item['output'])
                    counts.append(monitor(output, output.with_suffix('.log')))
                if counts:
                    diagnostics.update(native_runs=len(native), completed=sum(
                        item['resolved'] == 'finished' for item in native),
                                       done=sum(item[0] for item in counts),
                                       total=sum(item[1] for item in counts),
                                       estimated=True)
                    phase = ('departures' if any('dep_' in item['output'] for item in native)
                             else 'synthesis')
            except (OSError, KeyError, ValueError):
                pass
        return {'phase': phase, 'elapsed': max(
            0.0, (record.get('ended') or time.time()) - record['submitted']),
                'diagnostics': diagnostics}

    @property
    def logs(self):
        record = self._record()
        logs = list(record.get('logs') or ())
        if record.get('state') == 'running':
            try:
                from . import jobrecord
                logs.extend(str(Path(item['output']).with_suffix('.log'))
                            for item in jobrecord.scan(record['request']['output_root'])
                            if item.get('started', 0) >= record['submitted'] - 1
                            and Path(item['output']).with_suffix('.log').is_file())
            except (OSError, KeyError):
                pass
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
              'logs': [], 'result': None}
    _atomic(folder / RECORD, record)
    token = CancelToken()
    handle = WorkflowHandle(folder, token=token)

    def update(**values):
        current = _read(folder) or record
        current.update(values)
        _atomic(folder / RECORD, current)

    def work():
        update(state='running', phase='preparation', started=time.time())
        try:
            token.check()
            update(phase='synthesis')
            result = execute_request(resolved, cancel_token=token)
            update(phase='saving')
            result_path = result.save(folder / 'result')
            logs = [str(Path(run).with_suffix('.log')) for run in result.runs
                    if Path(run).with_suffix('.log').is_file()]
            update(state='finished', phase='finished', ended=time.time(),
                   result=str(result_path), runs=list(result.runs), logs=logs)
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

    future = WORKFLOWS.submit(work)
    handle._future = future
    with _lock:
        _live[identity] = handle
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
