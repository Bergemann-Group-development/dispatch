"""The calculations a durable workflow record can carry: the public operations.

A record names its operation as its task kind and carries the operation's JSON arguments
(``record['task']``); the driver that claims it looks the kind up here and runs it, so
every dispatch calculation a script or the GUI's tools submit (a synthesis, a fit, an
abundance fit, NLTE corrections, an ionisation balance, a departure solve) takes the same
detached path and survives the tab and the server. Each kind knows how to decode its
arguments, run them with a cancel token and a progress callback, save its result beside
the record and load it back, and how to encode a progress snapshot, so a document that
reopens the record while it runs can show the calculation where it is.

Arguments and results are the plain ``to_dict``/``from_dict`` forms the tumult classes
already provide for session snapshots; nothing here depends on Bokeh.
"""
from __future__ import annotations

import importlib
import json
from pathlib import Path

RESULT_DIR = 'result'
RESULT_FILE = 'result.json'       # its presence marks a record finished (submission)
PROGRESS_FILE = 'progress.json'


class TaskKind:
    """One kind of durable calculation. Subclasses override the four hooks; the default
    progress encoding is the result's own ``to_dict``."""
    name = ''
    progress = False        # whether run() reports intermediate results

    def decode(self, args):
        """JSON arguments -> the live objects ``run`` takes (validated)."""
        raise NotImplementedError

    def run(self, args, token, progress):
        raise NotImplementedError

    def save(self, result, folder):
        """Write ``result`` under ``folder`` (the record's ``result/``); return the
        manifest path."""
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / RESULT_FILE
        path.write_text(json.dumps({'kind': self.name, 'data': self.encode(result)},
                                   default=_json_default) + '\n')
        return path

    def load(self, folder):
        value = json.loads((Path(folder) / RESULT_FILE).read_text())
        return self.decode_result(value['data'])

    def encode(self, result):
        return result.to_dict()

    def decode_result(self, data):
        raise NotImplementedError


def _json_default(value):
    """numpy scalars and arrays, tuples of them, and anything else json refuses."""
    tolist = getattr(value, 'tolist', None)
    if tolist is not None:
        return tolist()
    item = getattr(value, 'item', None)
    if item is not None:
        return item()
    return str(value)


class Operation(TaskKind):
    """One public operation (``module.function``) on a
    :class:`~tumult.synthesis.SynthRequest`, run exactly as the scripting call would be.

    Arguments are ``{"request", "kwargs"}``: the request's ``to_dict`` and the operation's
    keyword arguments, with spectra, line-mask sets, atmospheres and measured results in
    their ``to_dict`` forms (a target or line mask may also be a file path)."""

    def __init__(self, name, module, function, progress, result_type):
        self.name, self.progress = name, progress
        self._module, self._function, self._result_type = module, function, result_type

    @property
    def function(self):
        return getattr(importlib.import_module(self._module), self._function)

    @staticmethod
    def encode_args(request, kwargs):
        values = {}
        for key, value in kwargs.items():
            values[key] = value.to_dict() if hasattr(value, 'to_dict') else \
                (str(value) if isinstance(value, Path) else value)
        return {'request': request.to_dict(), 'kwargs': values}

    def decode(self, args):
        from ..linemask import LineMaskSet
        from ..spectrum import Spectrum
        from .inputs import coerce_atmosphere
        from .request import SynthRequest
        kwargs = dict(args.get('kwargs') or {})
        if isinstance(kwargs.get('target'), dict):
            kwargs['target'] = Spectrum.from_dict(kwargs['target'])
        if isinstance(kwargs.get('line_mask'), dict):
            kwargs['line_mask'] = LineMaskSet.from_dict(kwargs['line_mask'])
        if isinstance(kwargs.get('reference_atmosphere'), dict):
            kwargs['reference_atmosphere'] = coerce_atmosphere(kwargs['reference_atmosphere'])
        if isinstance(kwargs.get('measured'), dict):
            kwargs['measured'] = self.decode_result(kwargs['measured'])
        return {'request': SynthRequest.from_dict(args['request']), 'kwargs': kwargs}

    def run(self, args, token, progress):
        extra = {'progress': progress} if self.progress else {}
        return self.function(args['request'], cancel_token=token, **extra, **args['kwargs'])

    def save(self, result, folder):
        if self._result_type == 'synthesis':
            return result.save(folder)
        return super().save(result, folder)

    def load(self, folder):
        if self._result_type == 'synthesis':
            from .request import SynthesisResult
            return SynthesisResult.load(Path(folder))
        return super().load(folder)

    def encode(self, result):
        if self._result_type == 'departures':
            from ..departures import to_record
            return to_record(result)
        return result.to_dict()

    def decode_result(self, data):
        if self._result_type == 'balance':
            from ..measurement import BalanceResult
            return BalanceResult.from_dict(data)
        if self._result_type == 'abundance':
            from ..abundance import AbundanceResult
            return AbundanceResult.from_dict(data)
        if self._result_type == 'fit':
            from .request import SpectrumFit
            return SpectrumFit.from_dict(data)
        if self._result_type == 'departures':
            from ..departures import from_record
            return from_record(data)
        raise ValueError(f'{self.name} has no progress or result record of this form')


OPERATIONS = (
    Operation('synthesis.synthesize', 'tumult.synthesis.operations', 'synthesize', False,
              'synthesis'),
    Operation('synthesis.fit', 'tumult.synthesis.operations', 'fit', False, 'fit'),
    Operation('departures.solve', 'tumult.departures', 'solve', False, 'departures'),
    Operation('abundance.fit', 'tumult.abundance', 'fit', True, 'abundance'),
    Operation('abundance.nlte_corrections', 'tumult.abundance', 'nlte_corrections', True,
              'abundance'),
    Operation('sparams.ionbalance', 'tumult.sparams', 'ionbalance', True, 'balance'),
)
KINDS = {kind.name: kind for kind in OPERATIONS}


def kind_of(record):
    """The :class:`TaskKind` a record carries. ``KeyError`` names an unknown kind, which
    includes records written before their operation had one."""
    name = (record.get('task') or {}).get('kind')
    try:
        return KINDS[name]
    except KeyError:
        raise KeyError(f'unknown workflow kind {name!r}') from None


def args_of(record):
    """The JSON arguments of a record's task."""
    return record['task']['args']


def write_progress(folder, kind, result):
    """Replace the record's progress snapshot with ``result``; a snapshot the reader
    cannot decode is worse than none, so a result that fails to encode is skipped."""
    path = Path(folder) / PROGRESS_FILE
    try:
        text = json.dumps({'kind': kind.name, 'data': kind.encode(result)},
                          default=_json_default)
    except (TypeError, ValueError):
        return
    tmp = path.with_name(f'{path.name}.partial')
    tmp.write_text(text + '\n')
    tmp.replace(path)


def read_progress(folder, kind):
    """The latest progress snapshot as a live result, or ``None`` (none written yet, or
    a half-written file being replaced)."""
    path = Path(folder) / PROGRESS_FILE
    try:
        value = json.loads(path.read_text())
        return kind.decode_result(value['data'])
    except (OSError, ValueError, KeyError, TypeError):
        return None
