"""Durable asynchronous submission of the same resolved calculation used synchronously."""
from __future__ import annotations

from copy import deepcopy
import json
import logging
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

    def __init__(self, folder, future=None, token=None, update=None):
        self.folder = Path(folder).expanduser().resolve()
        self.id = self.folder.name
        self._future = future
        self._token = token
        # writes terminal state for a workflow whose worker never runs (queued cancellation)
        self._update = update
        self._callbacks = []    # done callbacks owed after a queued cancellation

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
    def label(self):
        """The display name the submitter gave the workflow, else its id."""
        return self._record().get('label') or self.id

    @property
    def meta(self):
        """The submitter's display annotations (a JSON object), or ``None``."""
        return deepcopy(self._record().get('meta'))

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
        """Cancel this subscriber's calculation; never signal an unverified saved PID.

        A workflow still queued behind other pool work has its future cancelled, so the
        worker never starts; its record still settles as ``cancelled``. A running one is
        told through its token and settles when the worker notices."""
        if self._token is None or self._future is None or self._future.done():
            return False
        self._token.cancel()
        if self._future.cancel():
            self._settle_unstarted('cancelled')
        return True

    def _settle_unstarted(self, state, error=None):
        """Persist the terminal state of a workflow whose worker never ran, then tell
        its done callbacks; they observe the settled record, not the queued one."""
        values = dict(state=state, phase=state, ended=time.time())
        if error is not None:
            values['error'] = error
        with _lock:
            if self._update is not None:
                self._update(**values)
            _live.pop(self.id, None)
            owed, self._callbacks = self._callbacks, None
        for callback in owed:
            callback(self)

    @property
    def done(self):
        """Whether the worker has settled, or never will (opened without a live future)."""
        return self._future is None or self._future.done()

    def add_done_callback(self, callback):
        """Call ``callback(handle)`` once the workflow settles, from the worker thread
        or, for a queued cancellation, the cancelling one; immediately if it already has.
        Never blocks, and the record is terminal by the time the callback runs."""
        if self._future is None:
            callback(self)
            return
        with _lock:
            if self._callbacks is not None:
                self._callbacks.append(callback)
                # a worker-run future settles its own record before completing; a
                # cancelled one is reported by _settle_unstarted instead
                self._future.add_done_callback(
                    lambda future: None if future.cancelled() else callback(self))
                return
        callback(self)      # settled by a queued cancellation before registration

    def result(self, timeout=None):
        """The saved result. Raises :class:`CancelledError` for a cancelled workflow and
        ``RuntimeError`` for a failed or interrupted one; waits on a live future only if
        the workflow has not settled, so a done callback may call this freely."""
        if self._future is not None:
            try:
                self._future.result(timeout=timeout)
            except TimeoutError:
                raise
            except BaseException:
                pass    # cancelled or failed: the record below names the outcome
        state = self.status
        if state == 'cancelled':
            raise CancelledError(f'Workflow {self.id} was cancelled')
        if state != 'finished':
            raise RuntimeError(f'Workflow {self.id} is {state}: '
                               f'{self._record().get("error") or "no result"}')
        return SynthesisResult.load(self.folder / 'result')

    def to_dict(self):
        return self._record()


def submit(request, *, root=None, label=None, meta=None):
    """Copy, resolve and queue ``request``; return immediately with a durable handle.

    ``label`` and ``meta`` are the submitter's display name and annotations (JSON), kept
    in the record so a later process can list and present the workflow as it was
    submitted."""
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
              'result': None, 'label': label, 'meta': deepcopy(meta)}
    _atomic(folder / RECORD, record)
    token = CancelToken()

    def update(**values):
        with _lock:
            current = _read(folder) or record
            current.update(values)
            _atomic(folder / RECORD, current)

    handle = WorkflowHandle(folder, token=token, update=update)

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
            # the traceback belongs to the server log whichever document, if any, is
            # still waiting; the record keeps the concise cause for its Processes row
            logging.getLogger(__name__).exception('workflow %s failed', identity)
            update(state='failed', phase='failed', ended=time.time(),
                   error=f'{type(exc).__name__}: {exc}')
            raise
        finally:
            with _lock:
                _live.pop(identity, None)

    with _lock:
        _live[identity] = handle
        try:
            handle._future = WORKFLOWS.submit(work)
        except BaseException as exc:
            handle._settle_unstarted('failed', error=f'{type(exc).__name__}: {exc}')
            raise
    return handle


def cancel_live_workflows():
    """Cancel every workflow this process still holds live, for server shutdown; returns
    how many were told. Queued ones settle as ``cancelled`` at once, running ones when
    their worker notices the token."""
    with _lock:
        handles = list(_live.values())
    return sum(handle.cancel() for handle in handles)


def open_workflow(path, *, root=None):
    """Open a workflow by folder/path or stable id after a UI reconnect."""
    path = Path(path)
    folder = ((Path(root) / '_workflows' / path) if root is not None and len(path.parts) == 1
              else path).expanduser().resolve()
    with _lock:
        return _live.get(folder.name) or WorkflowHandle(folder)


def discover_workflows(root, limit=None):
    """Newest workflow handles under ``<root>/_workflows``, whatever their state: live
    ones (this process still holds their handle), finished results, and records left
    queued or running by a process that no longer exists, which report ``interrupted``.
    Only that folder is scanned, never the run cache beside it."""
    base = Path(root).expanduser().resolve() / '_workflows'
    folders = [path.parent for path in base.glob(f'*/{RECORD}')]
    folders.sort(key=lambda path: (_read(path) or {}).get('submitted', 0), reverse=True)
    handles = [open_workflow(path) for path in folders]
    return handles[:limit] if limit else handles


def run(request, *, cancel_token=None):
    """Blocking counterpart to :func:`submit`."""
    from .workflows import execute_request
    return execute_request(request, cancel_token=cancel_token)
