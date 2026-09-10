"""Execution invariants: shared budgets, cancellation, cache writers and processes."""
import concurrent.futures
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from tumult.synthesis.cancel import CancelToken, CancelledError
from tumult.synthesis.runtime import ExecutionService
from tumult.synthesis.executor import _run_batch


def test_global_budget_across_workflows():
    lock = threading.Lock()
    active = peak = 0
    barrier = threading.Event()

    def job(token, slots):
        nonlocal active, peak
        with lock:
            active += slots
            peak = max(peak, active)
        assert barrier.wait(3)
        with lock:
            active -= slots
        return slots

    with ExecutionService(3) as runtime:
        fs = [runtime.submit(str(i), '', i % 2, job, threads=2) for i in range(6)]
        barrier.set()
        assert all(f.result(3) == 2 for f in fs)
        assert peak <= 3


def test_shared_job_cancel_isolation():
    entered, release = threading.Event(), threading.Event()
    calls = []

    def job(token, slots):
        calls.append(1)
        entered.set()
        assert release.wait(3)
        token.check()
        return 42

    with ExecutionService(2) as runtime:
        a, b = CancelToken(), CancelToken()
        f1 = runtime.submit('same', '', 'a', job, token=a)
        assert entered.wait(3)
        f2 = runtime.submit('same', '', 'b', job, token=b)
        a.cancel()
        with pytest.raises(CancelledError):
            f1.result(1)
        release.set()
        assert f2.result(3) == 42
        assert len(calls) == 1


def test_cancel_queued_and_running_jobs_leaves_service_usable():
    entered = threading.Event()
    token = CancelToken()

    def job(job_token, slots):
        entered.set()
        stopped = threading.Event()
        detach = job_token.add_callback(stopped.set)
        try:
            assert stopped.wait(3)
            job_token.check()
        finally:
            detach()

    with ExecutionService(1) as runtime:
        running = runtime.submit('a', '', 'a', job, token=token)
        assert entered.wait(3)
        queued = runtime.submit('b', '', 'b', lambda *_: pytest.fail('cancelled job ran'), token=token)
        token.cancel()
        for f in (running, queued):
            with pytest.raises(CancelledError):
                f.result(2)
        assert runtime.submit('c', '', 'c', lambda *_: 7).result(3) == 7


def test_exception_releases_capacity_and_conflicting_output_rejected():
    entered, release = threading.Event(), threading.Event()

    def fail(*_):
        entered.set()
        assert release.wait(3)
        raise ValueError('numerical failure')

    with ExecutionService(1) as runtime:
        f = runtime.submit('same', 'first', 'a', fail)
        assert entered.wait(3)
        with pytest.raises(ValueError, match='disagree'):
            runtime.submit('same', 'different', 'b', fail)
        release.set()
        with pytest.raises(ValueError, match='numerical'):
            f.result(3)
        assert runtime.submit('next', '', 'b', lambda *_: 1).result(3) == 1


def test_batch_deduplicates_and_preserves_keys(monkeypatch, tmp_path):
    from tumult.synthesis import executor
    entered, release = threading.Event(), threading.Event()
    calls = []

    def run(key, name, nml, kwargs):
        calls.append(name)
        entered.set()
        assert release.wait(3)
        return key, 'spectrum'

    monkeypatch.setattr(executor, '_run_single', run)
    jobs = [(1, tmp_path / 'same', {}), (2, tmp_path / 'same', {})]
    with ExecutionService(2) as runtime, concurrent.futures.ThreadPoolExecutor(1) as caller:
        f = caller.submit(_run_batch, jobs, {}, 2, runtime=runtime)
        assert entered.wait(3)
        # Both submissions are made before the bounded batch starts waiting.
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline:
            with runtime._lock:
                if len(runtime._tasks[str((tmp_path / 'same').resolve())].subscribers) == 2:
                    break
            time.sleep(.005)
        release.set()
        assert f.result(3) == {1: 'spectrum', 2: 'spectrum'}
        assert len(calls) == 1


def test_launch_cancellation_kills_child_and_invalidates_partial_cache(tmp_path):
    from tumult import launch
    output = tmp_path / 'job'
    ready = tmp_path / 'ready'
    executable = tmp_path / 'dispatch.x'
    executable.write_text(
        f'#!{sys.executable}\n'
        'from pathlib import Path\nimport time\n'
        f'p=Path({str(output)!r}); p.mkdir(parents=True, exist_ok=True)\n'
        '(p / "abund.txt").write_text("partial")\n'
        f'Path({str(ready)!r}).touch()\n'
        'time.sleep(30)\n')
    executable.chmod(0o755)
    token = CancelToken()
    with concurrent.futures.ThreadPoolExecutor(1) as executor:
        f = executor.submit(launch, str(tmp_path / 'job'),
                            {'io_params': {'datadir': str(tmp_path / 'data')}},
                            cwd=tmp_path, verbose=0, cancel_token=token)
        deadline = time.monotonic() + 5
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(.01)
        assert ready.exists()
        token.cancel()
        with pytest.raises(CancelledError):
            f.result(3)
        assert not (output / 'abund.txt').exists()


def test_round_robin_between_workflows():
    entered, release = threading.Event(), threading.Event()
    order = []

    def blocker(*_):
        entered.set()
        assert release.wait(3)

    def job(name):
        def run(*_):
            order.append(name)
        return run

    with ExecutionService(1) as runtime:
        first = runtime.submit('block', '', 'a', blocker)
        assert entered.wait(3)
        fs = [runtime.submit('a1', '', 'a', job('a1')),
              runtime.submit('a2', '', 'a', job('a2')),
              runtime.submit('b1', '', 'b', job('b1'))]
        release.set()
        first.result(3)
        for f in fs:
            f.result(3)
        assert order == ['a1', 'b1', 'a2']


def test_shutdown_cancels_running_and_queued():
    entered = threading.Event()
    runtime = ExecutionService(1)

    def run(token, _):
        stop = threading.Event()
        detach = token.add_callback(stop.set)
        entered.set()
        try:
            assert stop.wait(3)
            token.check()
        finally:
            detach()

    a = runtime.submit('a', '', 'a', run)
    assert entered.wait(3)
    b = runtime.submit('b', '', 'b', lambda *_: pytest.fail('queued job ran'))
    runtime.shutdown()
    for f in (a, b):
        with pytest.raises(CancelledError):
            f.result(1)
    with pytest.raises(RuntimeError, match='shut down'):
        runtime.submit('c', '', 'c', run)


def test_explicit_scope_used_by_public_compute_calls():
    from tumult.synthesis import execution_scope
    from tumult.synthesis.runtime import get_runtime
    default = get_runtime()
    with execution_scope(2) as runtime:
        assert get_runtime() is runtime
        assert runtime.budget == 2
    assert get_runtime() is default


def test_concurrent_batches_share_real_subprocess_and_cancel_independently(tmp_path, monkeypatch):
    import tumult as tm
    executable = tmp_path / 'dispatch.x'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import sys, os, time\nfrom pathlib import Path\n'
        'p=Path(sys.argv[1]).with_suffix(""); p.mkdir(parents=True, exist_ok=True)\n'
        'with (p / "starts").open("a") as f: f.write(os.environ["OMP_NUM_THREADS"] + "\\n")\n'
        'time.sleep(.3)\n'
        '(p / "abund.txt").write_text("finished")\n')
    executable.chmod(0o755)
    monkeypatch.setattr(tm, 'read', lambda path: path.name)
    jobs = [(i, tmp_path / f'job{i}', {}) for i in range(3)]
    a, b = CancelToken(), CancelToken()
    with ExecutionService(2) as runtime, concurrent.futures.ThreadPoolExecutor(2) as callers:
        fa = callers.submit(_run_batch, jobs, {'cwd': tmp_path, 'verbose': 0}, 2, a, runtime=runtime)
        fb = callers.submit(_run_batch, jobs, {'cwd': tmp_path, 'verbose': 0}, 2, b, runtime=runtime)
        deadline = time.monotonic() + 3
        shared = False
        while time.monotonic() < deadline:
            with runtime._lock:
                shared = any(len(t.subscribers) == 2 for t in runtime._tasks.values())
            if shared:
                break
            time.sleep(.005)
        assert shared
        a.cancel()
        with pytest.raises(CancelledError):
            fa.result(3)
        assert fb.result(5) == {i: f'job{i}' for i in range(3)}
        for i in range(3):
            assert (tmp_path / f'job{i}' / 'starts').read_text() == '1\n'


def test_compute_spectra_folds_abundance_only_batch(monkeypatch, tmp_path):
    from tumult.synthesis import api
    submitted = []

    class FakeRun:
        has_spectrum = True
        mode = 'LTE'
        folder = tmp_path / 'run'

        def crop_set(self, index, *args, **kwargs):
            return index

    monkeypatch.setattr(api, 'select_nodes', lambda *a, **k: (
        [('node', 5000., 4., 0., 1., 1., 'atmos', 'absmet')], 1.))

    def run_batch(jobs, *args, **kwargs):
        submitted.append(list(jobs))
        return {key: FakeRun() for key, _path, _nml in jobs}

    monkeypatch.setattr(api, '_run_batch', run_batch)
    monkeypatch.setattr(api, '_interp_crop',
                        lambda grid, *args: grid.spectra[0].index)
    settings = [
        {'composition_params': {'abundances': 'Mg=-0.2 Si=-0.2',
                                'relative_abundances': True}},
        {'composition_params': {'abundances': 'Mg=+0.2 Si=+0.2',
                                'relative_abundances': True}},
    ]
    result = api.compute_spectra(
        tmp_path, (5000., 4., 0., 1.), {}, settings_list=settings,
        max_threads=1, xmin=5000., xmax=5001., nlte=False)
    assert result == [2, 3]
    assert [len(batch) for batch in submitted] == [1]
    nml = submitted[0][0][2]
    manifest = Path(nml['output_params']['set_file'])
    assert manifest.is_file()
    assert 'Mg' in manifest.read_text() and 'Si' in manifest.read_text()


def test_compute_spectra_keeps_n_run_fallback(monkeypatch, tmp_path):
    from tumult.synthesis import api
    submitted = []

    class FakeRun:
        has_spectrum = True
        mode = 'LTE'
        folder = tmp_path / 'run'

    monkeypatch.setattr(api, 'select_nodes', lambda *a, **k: (
        [('node', 5000., 4., 0., 1., 1., 'atmos', 'absmet')], 1.))

    def run_batch(jobs, *args, **kwargs):
        submitted.append(list(jobs))
        return {key: FakeRun() for key, _path, _nml in jobs}

    monkeypatch.setattr(api, '_run_batch', run_batch)
    monkeypatch.setattr(api, '_interp_crop', lambda *args: 'spectrum')
    settings = [
        {'composition_params': {'abundances': 'Fe=-0.2'},
         'linelist_params': {'only_elements': 'Fe'}},
        {'composition_params': {'abundances': 'Fe=+0.2'},
         'linelist_params': {'only_elements': 'Mg'}},
    ]
    result = api.compute_spectra(
        tmp_path, (5000., 4., 0., 1.), {}, settings_list=settings,
        max_threads=2, xmin=5000., xmax=5001., nlte=False)
    assert result == ['spectrum', 'spectrum']
    assert [len(batch) for batch in submitted] == [2]


def test_provenance_round_trips_through_the_cache_listing(tmp_path, monkeypatch):
    """A run states its own point; the listing reads it back rather than re-deriving it."""
    from tumult.synthesis import api
    from tumult.synthesis import workflows

    run_dir = tmp_path / 'data' / 'synth' / 'w5000_5010_abcd1234' / 't58g44m+00v10'
    run_dir.mkdir(parents=True)
    (run_dir / 'params.nml').write_text("&atmos_params /\n")
    runs = {'node': SimpleNamespace(folder=run_dir)}
    nodes = [('node', 5750.0, 4.4, 0.0, 1.0, 1.0, 'atmos', 'absmet')]
    api._save_provenance(runs, nodes, {'composition_params': {
        'abundances': 'Fe=-0.20', 'relative_abundances': True}})

    monkeypatch.chdir(tmp_path)
    entry, = workflows.list_cache({'synth_output_dir': 'synth'})
    # the folder name rounds 5750 K to 't58'; the record does not
    assert entry['star'] == {'teff': 5750.0, 'logg': 4.4, 'feh': 0.0, 'vmic': 1.0}
    assert entry['abund'] == {'items': [('Fe', -0.2)], 'relative': True}
    assert entry['window'] == (5000, 5010)


def test_cache_listing_falls_back_to_the_folder_name(tmp_path, monkeypatch):
    """A run from before provenance.json keeps listing, at the folder name's resolution."""
    from tumult.synthesis import workflows

    run_dir = tmp_path / 'data' / 'synth' / 'w5000_5010_abcd1234' / 't58g44m+00v10'
    run_dir.mkdir(parents=True)
    (run_dir / 'params.nml').write_text("&atmos_params /\n")

    monkeypatch.chdir(tmp_path)
    entry, = workflows.list_cache({'synth_output_dir': 'synth'})
    assert entry['star'] == {'teff': 58, 'logg': 4.4, 'feh': 0.0, 'vmic': 1.0}
    assert entry['abund'] == {'items': [], 'relative': False}


def test_launch_runs_relative_inputs_in_execution_directory(tmp_path, monkeypatch):
    from tumult import launch
    execution, client = tmp_path/'execution', tmp_path/'client'
    execution.mkdir(); client.mkdir()
    executable = execution/'dispatch.x'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import sys\nfrom pathlib import Path\n'
        'assert Path(sys.argv[1]).is_file()\n'
        'out=Path("data")/Path(sys.argv[1]).with_suffix("")\n'
        'out.mkdir(parents=True)\n(out/"abund.txt").write_text("done")\n')
    executable.chmod(0o755)
    monkeypatch.chdir(client)
    output, ran = launch('cache/job', cwd=execution, verbose=0)
    assert ran and output == execution/'data/cache/job'
    assert (execution/'cache/job.nml').exists()
    assert not (client/'cache').exists()
    assert launch('cache/job', cwd=execution, verbose=0) == (output, False)
