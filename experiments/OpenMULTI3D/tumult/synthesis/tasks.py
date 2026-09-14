"""The kinds of calculation a durable workflow record can carry.

A record names its kind and JSON arguments (``record['task']``); the driver that claims
it looks the kind up here and runs it, so every dispatch calculation the GUI's tools
start (a synthesis, a fit, an abundance run, an ionisation balance, NLTE corrections)
takes the same detached path and survives the tab and the server. Each kind knows how
to decode its arguments, run them with a cancel token and a progress callback, save its
result beside the record and load it back, and how to encode a progress snapshot, so a
document that reopens the record while it runs can show the calculation where it is.

Arguments and results are the plain ``to_dict``/``from_dict`` forms the tumult classes
already provide for session snapshots; nothing here depends on Bokeh.
"""
from __future__ import annotations

from dataclasses import asdict
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


# --- synthesis: the original workflow record --------------------------------------------
class Synthesis(TaskKind):
    """One resolved synthesis request; the result is the spectra it produced."""
    name = 'synthesis'

    def decode(self, args):
        from .request import ResolvedSynthRequest
        return ResolvedSynthRequest.from_dict(args)

    def run(self, resolved, token, progress):
        from .workflows import execute_request
        return execute_request(resolved, cancel_token=token)

    def save(self, result, folder):
        return result.save(folder)

    def load(self, folder):
        from .request import SynthesisResult
        return SynthesisResult.load(Path(folder))


# --- the tool workflows ---------------------------------------------------------------------
class NLTECorrections(TaskKind):
    """Stage 1 (the 1D LTE reference and its measured lines) and stage 2 (the NLTE curve
    of growth inverted to each line's correction) of :mod:`tumult.abundance`. With a
    stage-1 ``result`` in the arguments only stage 2 runs. The stage-1 result is reported
    as progress, so a document sees the reference lines before the corrections."""
    name = 'nlte_corrections'
    progress = True

    def decode(self, args):
        from ..abundance import AbundanceResult, CorrectionOpts
        from .request import NativeRequest
        return {'request': NativeRequest.from_dict(args['request']),
                'atom': str(args['atom']),
                'ref_grid': dict(args['ref_grid'] or {}),
                'opts': CorrectionOpts(**dict(args.get('opts') or {})),
                'result': (AbundanceResult.from_dict(args['result'])
                           if args.get('result') else None)}

    def run(self, args, token, progress):
        from ..abundance import measure_reference_lines, run_nlte_corrections
        request, atom, opts = args['request'], args['atom'], args['opts']
        result = args['result']
        if result is None:
            result = measure_reference_lines(request, atom, ref_grid=args['ref_grid'],
                                             opts=opts, cancel_token=token)
            result.synthesis = asdict(request)      # what names it; stage 2 rewrites it
            progress(result)
            token.check()           # a cancel during the reference must not start stage 2
        return run_nlte_corrections(request, result, atom, opts=opts, progress=progress,
                                    cancel_token=token)

    def decode_result(self, data):
        from ..abundance import AbundanceResult
        return AbundanceResult.from_dict(data)


class Abundance(TaskKind):
    """The per-line synthesis and fit of a measured line set (:func:`run_abundance`),
    then the error budget (:func:`finalize_result`)."""
    name = 'abundance'
    progress = True

    def decode(self, args):
        from ..abundance import AbundanceResult
        from ..spectrum import Spectrum
        from .request import NativeRequest
        return {'request': NativeRequest.from_dict(args['request']),
                'result': AbundanceResult.from_dict(args['result']),
                'observation': Spectrum.from_dict(args['observation']),
                'snr': None if args.get('snr') is None else float(args['snr']),
                'cont_pct': float(args.get('cont_pct', float('nan')))}

    def run(self, args, token, progress):
        from ..abundance import finalize_result, run_abundance
        result, obs = args['result'], args['observation']
        run_abundance(args['request'], result, obs, cancel_token=token, progress=progress)
        finalize_result(result, obs, args['snr'], args['cont_pct'])
        return result

    def decode_result(self, data):
        from ..abundance import AbundanceResult
        return AbundanceResult.from_dict(data)


class IonisationBalance(TaskKind):
    """Measure a mask set's windows on the observation (when ``masks`` are given) and
    iterate the excitation-ionisation balance (:mod:`tumult.sparams`)."""
    name = 'ionisation_balance'
    progress = True

    def decode(self, args):
        from ..linemask import LineMask
        from ..measurement import BalanceResult
        from ..spectrum import Spectrum
        from .request import NativeRequest
        masks = args.get('masks')
        return {'request': NativeRequest.from_dict(args['request']),
                'result': BalanceResult.from_dict(args['result']),
                'masks': None if masks is None else [LineMask.from_dict(m) for m in masks],
                'observation': (Spectrum.from_dict(args['observation'])
                                if args.get('observation') else None),
                'linelist': args.get('linelist'),
                'opts': dict(args.get('opts') or {})}

    def run(self, args, token, progress):
        from .. import sparams
        result = args['result']
        if args['masks'] is not None:
            sparams.measure_balance_lines(result, args['masks'], args['observation'],
                                          args['linelist'])
            progress(result)
        return sparams.solve_ionisation_balance(args['request'], result, args['opts'],
                                                cancel_token=token, progress=progress)

    def decode_result(self, data):
        from ..measurement import BalanceResult
        return BalanceResult.from_dict(data)


class Fit(TaskKind):
    """An observation fitted by adaptive synthesis (:func:`workflows.fit`); the result
    is ``(best, raw, post)``: the fitted values, the best-fit spectrum without
    post-processing and the broadening to apply as a live convolution."""
    name = 'fit'

    def decode(self, args):
        from ..spectrum import Spectrum
        from .request import FitRequest
        return {'observation': Spectrum.from_dict(args['observation']),
                'request': FitRequest.from_dict(args['request'])}

    def run(self, args, token, progress):
        from .workflows import fit
        return fit(args['observation'], args['request'], cancel_token=token)

    def encode(self, result):
        best, raw, post = result
        return {'best': {str(k): float(v) for k, v in dict(best).items()},
                'raw': raw.to_dict(), 'post': {str(k): float(v) for k, v in dict(post).items()}}

    def decode_result(self, data):
        from ..spectrum import Spectrum
        return dict(data['best']), Spectrum.from_dict(data['raw']), dict(data['post'])


# --- the public operations ---------------------------------------------------------------
OPERATIONS = {
    'synthesis.synthesize': ('tumult.synthesis.operations', 'synthesize', False),
    'departures.solve': ('tumult.departures', 'solve', False),
    'abundance.fit': ('tumult.abundance', 'fit', True),
    'abundance.nlte_corrections': ('tumult.abundance', 'nlte_corrections', True),
    'sparams.ionbalance': ('tumult.sparams', 'ionbalance', True),
}


def _operation(name):
    import importlib
    try:
        module, function, progress = OPERATIONS[name]
    except KeyError:
        raise ValueError(f'unknown operation {name!r}; one of {sorted(OPERATIONS)}') from None
    return getattr(importlib.import_module(module), function), progress


class Operation(TaskKind):
    """One public scripting operation (:data:`OPERATIONS`) on a
    :class:`~tumult.synthesis.SynthRequest`, run exactly as the scripting call would be.

    Arguments are ``{"operation", "request", "kwargs"}``: the request's ``to_dict`` and the
    operation's keyword arguments, with a target spectrum and a line-mask set in their
    ``to_dict`` forms (or as file paths)."""
    name = 'operation'
    progress = True

    @staticmethod
    def encode_args(operation, request, kwargs):
        _operation(operation)
        values = {}
        for key, value in kwargs.items():
            values[key] = value.to_dict() if hasattr(value, 'to_dict') else \
                (str(value) if isinstance(value, Path) else value)
        return {'operation': operation, 'request': request.to_dict(), 'kwargs': values}

    def decode(self, args):
        from ..linemask import LineMaskSet
        from ..spectrum import Spectrum
        from .request import SynthRequest
        function, progress = _operation(args['operation'])
        kwargs = dict(args.get('kwargs') or {})
        if isinstance(kwargs.get('target'), dict):
            kwargs['target'] = Spectrum.from_dict(kwargs['target'])
        if isinstance(kwargs.get('line_mask'), dict):
            kwargs['line_mask'] = LineMaskSet.from_dict(kwargs['line_mask'])
        return {'operation': args['operation'], 'function': function, 'progress': progress,
                'request': SynthRequest.from_dict(args['request']), 'kwargs': kwargs}

    def run(self, args, token, progress):
        extra = {'progress': progress} if args['progress'] else {}
        return args['function'](args['request'], cancel_token=token, **extra,
                                **args['kwargs'])

    def save(self, result, folder):
        from .request import SynthesisResult
        if isinstance(result, SynthesisResult):
            return result.save(folder)
        return super().save(result, folder)

    def load(self, folder):
        from .request import SynthesisResult
        value = json.loads((Path(folder) / RESULT_FILE).read_text())
        if 'native_request' in value:
            return SynthesisResult.load(Path(folder))
        return self.decode_result(value['data'])

    def encode(self, result):
        from ..departures import to_record
        from ..measurement import BalanceResult
        if isinstance(result, BalanceResult):
            return {'type': 'balance', 'data': result.to_dict()}
        if hasattr(result, 'to_abundances'):
            return {'type': 'abundance', 'data': result.to_dict()}
        return {'type': 'departures', 'data': to_record(result)}

    def decode_result(self, data):
        from ..abundance import AbundanceResult
        from ..departures import from_record
        from ..measurement import BalanceResult
        kind = data['type']
        if kind == 'balance':
            return BalanceResult.from_dict(data['data'])
        if kind == 'abundance':
            return AbundanceResult.from_dict(data['data'])
        return from_record(data['data'])


KINDS = {kind.name: kind for kind in (Synthesis(), NLTECorrections(), Abundance(),
                                      IonisationBalance(), Fit(), Operation())}


def kind_of(record):
    """The :class:`TaskKind` a record carries; a record without a ``task`` block is a
    synthesis (the original record layout). ``KeyError`` names an unknown kind."""
    name = ((record.get('task') or {}).get('kind')) or Synthesis.name
    try:
        return KINDS[name]
    except KeyError:
        raise KeyError(f'unknown workflow kind {name!r}') from None


def args_of(record):
    """The JSON arguments of a record's task: ``record['request']`` for a synthesis."""
    task = record.get('task') or {}
    if task.get('kind', Synthesis.name) == Synthesis.name:
        return record['request']
    return task['args']


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
