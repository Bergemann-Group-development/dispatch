"""Execution invariants: shared budgets, cancellation, cache writers and processes."""
import concurrent.futures
import os
import socket
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


def test_balanced_allocation_gives_a_single_ready_job_the_whole_budget(tmp_path, monkeypatch):
    """plan/multi-processing.md's ready-queue table: one ready job on ``F`` free workers
    gets all ``F`` of them, not a fixed per-job thread count."""
    import tumult as tm
    executable = tmp_path / 'dispatch.x'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import sys, os\nfrom pathlib import Path\n'
        'p=Path(sys.argv[1]).with_suffix(""); p.mkdir(parents=True, exist_ok=True)\n'
        '(p / "starts").write_text(os.environ["OMP_NUM_THREADS"])\n'
        '(p / "abund.txt").write_text("finished")\n')
    executable.chmod(0o755)
    monkeypatch.setattr(tm, 'read', lambda path: path.name)
    jobs = [(0, tmp_path / 'job0', {})]
    with ExecutionService(4) as runtime:
        result = _run_batch(jobs, {'cwd': tmp_path, 'verbose': 0}, 4, runtime=runtime)
    assert result == {0: 'job0'}
    assert (tmp_path / 'job0' / 'starts').read_text() == '4'


def test_balanced_allocation_splits_free_budget_across_ready_jobs(tmp_path, monkeypatch):
    """plan/multi-processing.md's ready-queue table: two jobs ready on seven free workers
    split ``4, 3`` (``base = 7 // 2 = 3``, the first job getting the one-slot remainder)."""
    import tumult as tm
    executable = tmp_path / 'dispatch.x'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import sys, os\nfrom pathlib import Path\n'
        'p=Path(sys.argv[1]).with_suffix(""); p.mkdir(parents=True, exist_ok=True)\n'
        '(p / "starts").write_text(os.environ["OMP_NUM_THREADS"])\n'
        '(p / "abund.txt").write_text("finished")\n')
    executable.chmod(0o755)
    monkeypatch.setattr(tm, 'read', lambda path: path.name)
    jobs = [(i, tmp_path / f'job{i}', {}) for i in range(2)]
    with ExecutionService(7) as runtime:
        result = _run_batch(jobs, {'cwd': tmp_path, 'verbose': 0}, 7, runtime=runtime)
    assert result == {0: 'job0', 1: 'job1'}
    assert (tmp_path / 'job0' / 'starts').read_text() == '4'
    assert (tmp_path / 'job1' / 'starts').read_text() == '3'


def test_balanced_allocation_grants_freed_slots_to_a_job_unlocked_later(monkeypatch, tmp_path):
    """With more ready jobs than free workers, the extra ones wait; a completion frees its
    slots and the next ready job is admitted on the very same accounting, never exceeding
    the budget (plan/multi-processing.md: "On every completion, release that job's slots,
    update the dependency graph and immediately allocate the free slots")."""
    from tumult.synthesis import executor
    entered = {name: threading.Event() for name in ('job0', 'job1', 'job2', 'job3')}
    release = {name: threading.Event() for name in ('job0', 'job1', 'job2', 'job3')}
    seen = {}

    def run(key, name, nml, kwargs):
        name = str(name.name)
        seen[name] = kwargs['threads']
        entered[name].set()
        assert release[name].wait(3)
        return key, name

    monkeypatch.setattr(executor, '_run_single', run)
    jobs = [(i, tmp_path / f'job{i}', {}) for i in range(4)]
    with ExecutionService(3) as runtime, concurrent.futures.ThreadPoolExecutor(1) as caller:
        f = caller.submit(_run_batch, jobs, {}, 3, runtime=runtime)
        for name in ('job0', 'job1', 'job2'):
            assert entered[name].wait(3)
        assert not entered['job3'].is_set()             # only 3 workers: it waits
        assert seen == {'job0': 1, 'job1': 1, 'job2': 1}   # n=3 <= F=3: one thread each
        release['job0'].set()
        assert entered['job3'].wait(3)                   # job0's freed slot admits it
        assert seen['job3'] == 1
        for name in ('job1', 'job2', 'job3'):
            release[name].set()
        assert f.result(3) == {0: 'job0', 1: 'job1', 2: 'job2', 3: 'job3'}


def test_cache_hit_gets_no_allocation_and_does_not_shrink_others(monkeypatch, tmp_path):
    """plan/multi-processing.md: "Identify valid disk-cache hits ... load them without
    assigning native CPU slots" and cache hits must not reduce uncached jobs' share."""
    from tumult.synthesis import executor
    cached = tmp_path / 'cached'
    (cached / 'inner').mkdir(parents=True)
    (cached / 'inner' / 'abund.txt').write_text('done')     # already a finished run
    entered = threading.Event()
    seen = {}

    def run(key, name, nml, kwargs):
        if name.name == 'inner':
            return key, 'cached-result'
        seen['fresh'] = kwargs['threads']
        entered.set()
        return key, 'fresh-result'

    monkeypatch.setattr(executor, '_run_single', run)
    jobs = [(0, cached / 'inner', {}), (1, tmp_path / 'fresh', {})]
    with ExecutionService(5) as runtime:
        result = _run_batch(jobs, {'cwd': tmp_path}, 5, runtime=runtime)
    assert result == {0: 'cached-result', 1: 'fresh-result'}
    assert entered.is_set()
    assert seen['fresh'] == 5      # the sole uncached job gets the full free budget


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
        # exist_ok, like dispatch.x's own mkdir -p: the launcher writes the run's
        # job record into this folder before starting the executable.
        'out.mkdir(parents=True, exist_ok=True)\n(out/"abund.txt").write_text("done")\n')
    executable.chmod(0o755)
    monkeypatch.chdir(client)
    output, ran = launch('cache/job', cwd=execution, verbose=0)
    assert ran and output == execution/'data/cache/job'
    assert (execution/'cache/job.nml').exists()
    assert not (client/'cache').exists()
    assert launch('cache/job', cwd=execution, verbose=0) == (output, False)


def _stub_dispatch(execution, body):
    """A fake dispatch.x that runs ``body`` (Python source) in the execution directory."""
    executable = execution / 'dispatch.x'
    executable.write_text(f'#!{sys.executable}\n'
                          'import sys, time\nfrom pathlib import Path\n'
                          'out=Path("data")/Path(sys.argv[1]).with_suffix("")\n'
                          'out.mkdir(parents=True, exist_ok=True)\n' + body)
    executable.chmod(0o755)
    return executable


def test_job_record_states_a_finished_run(tmp_path, monkeypatch):
    """A local run records itself at launch and resolves as finished afterwards."""
    from tumult import launch
    from tumult.synthesis import jobrecord

    execution = tmp_path / 'execution'
    execution.mkdir()
    _stub_dispatch(execution, '(out/"abund.txt").write_text("done")\n')
    monkeypatch.chdir(execution)

    output, _ = launch('cache/job', {'atmos_params': {'nx': 7}}, verbose=0)
    record = jobrecord.read(output)
    assert record['state'] == 'finished' and record['pid']
    assert record['settings'] == jobrecord.settings_hash({'atmos_params': {'nx': 7}})
    assert jobrecord.state_of(output) == 'finished'
    listed, = jobrecord.scan(execution / 'data')
    assert listed['resolved'] == 'finished' and Path(listed['output']) == output


def test_job_record_survives_a_crashed_run(tmp_path, monkeypatch):
    """A run that never wrote abund.txt is reported failed, not silently absent."""
    from tumult import launch
    from tumult.synthesis import jobrecord

    execution = tmp_path / 'execution'
    execution.mkdir()
    _stub_dispatch(execution, 'sys.exit(3)\n')
    monkeypatch.chdir(execution)

    output, _ = launch('cache/job', verbose=0)
    assert jobrecord.state_of(output) == 'failed'
    # and the failure is not mistaken for a cache hit on the next launch
    assert launch('cache/job', verbose=0)[1] is True


def test_job_record_reports_a_killed_run_as_lost(tmp_path):
    """A record left behind by a server restart resolves from the process table."""
    from tumult.synthesis import jobrecord

    folder = tmp_path / 'run'
    jobrecord.write(folder, {'atmos_params': {}}, pid=2 ** 30)   # no such process
    assert jobrecord.state_of(folder) == 'lost'
    jobrecord.write(folder, {'atmos_params': {}}, pid=os.getpid())
    assert jobrecord.state_of(folder) == 'running'
    (folder / 'abund.txt').write_text('done')
    assert jobrecord.state_of(folder) == 'finished'  # the cache outranks the record


def _native_record(folder, **fields):
    """What SRC/job_record_mod.f90 writes, with the given fields overriding."""
    values = {'state': 'running', 'pid': os.getpid(), 'host': socket.gethostname(),
              'cwd': '/work', 'input': 'run.nml', 'started': '2026-09-13T11:15:20+02:00',
              'ended': '', 'mpi_size': 1, 'omp_threads': 4, 'slurm_job_id': '',
              'dispatch_hash': 'c098fe60', **fields}
    lines = ['&job_record']
    for key, value in values.items():
        lines.append(f' {key}={value!r},' if isinstance(value, str) else f' {key}={value},')
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'job.nml').write_text('\n'.join(lines) + '\n /\n')


def test_dispatch_own_record_makes_a_shell_run_discoverable(tmp_path):
    """A run started outside the launcher has only job.nml; it is listed, attributed
    and resolved from the process table like a launched one."""
    from tumult.synthesis import jobrecord
    folder = tmp_path / 'data' / 'dep_x'
    _native_record(folder)
    record = jobrecord.read(folder)
    assert record['pid'] == os.getpid() and record['input'] == 'run.nml'
    assert record['dispatch_hash'] == 'c098fe60' and record['slurm_id'] is None
    assert jobrecord.state_of(folder) == 'running'
    listed, = jobrecord.scan(tmp_path / 'data')
    assert listed['resolved'] == 'running' and Path(listed['output']) == folder
    _native_record(folder, pid=2 ** 30)
    assert jobrecord.state_of(folder) == 'lost'
    _native_record(folder, pid=2 ** 30, host='elsewhere')
    assert jobrecord.state_of(folder) == 'running'      # not inspectable from here
    _native_record(folder, state='finished', ended='2026-09-13T11:15:29+02:00')
    assert jobrecord.state_of(folder) == 'finished'     # a normal exit, even before abund.txt


def test_launcher_and_dispatch_records_merge(tmp_path, monkeypatch):
    """The launcher's terminal states outrank dispatch's stale 'running'; dispatch's pid
    outranks the launcher's, which under sbatch or mpirun is not the run's."""
    from tumult.synthesis import jobrecord, scheduler
    folder = tmp_path / 'run'
    jobrecord.write(folder, {'atmos_params': {}}, pid=12345)
    _native_record(folder)
    assert jobrecord.read(folder)['pid'] == os.getpid()
    assert jobrecord.read(folder)['launcher_pid'] == 12345
    jobrecord.finish(folder, 'cancelled')
    assert jobrecord.state_of(folder) == 'cancelled'
    # a Slurm run is resolved through the scheduler, whichever writer named the job
    _native_record(folder, slurm_job_id='4242')
    jobrecord.write(folder, {'atmos_params': {}}, state='submitted')
    monkeypatch.setattr(scheduler, 'state', lambda job_id: 'running' if job_id == '4242' else None)
    assert jobrecord.state_of(folder) == 'running'
    monkeypatch.setattr(scheduler, 'state', lambda job_id: None)
    assert jobrecord.state_of(folder) == 'lost'


def test_job_record_marks_a_cancelled_run(tmp_path, monkeypatch):
    """Cancelling a run leaves a record saying so, not one that still claims to run."""
    from tumult import launch
    from tumult.synthesis import jobrecord

    execution = tmp_path / 'execution'
    execution.mkdir()
    _stub_dispatch(execution, 'time.sleep(30)\n')
    monkeypatch.chdir(execution)

    token = CancelToken()
    threading.Timer(0.4, token.cancel).start()
    with pytest.raises(CancelledError):
        launch('cache/job', verbose=0, cancel_token=token)
    assert jobrecord.state_of(execution / 'data/cache/job') == 'cancelled'


MONITOR_LOG = """\
 ============================== LTE Iteration  ==================================
  local =  2.7826-320, global =  0.0000E+00, lvl =   0, pos=[   0   0   0] Mg
 ============================== Iteration  1/99 =================================
  local =  1.2475E-04, global = -2.3649E-01, lvl =  86, pos=[   2   1 128] Mg
  flux global =  5.8673E-01, i_nu =   1, lambda =   3112.49 AA
 ============================== Iteration  2/99 =================================
  local =  8.8159E-05, global = -9.1297E-02, lvl =  86, pos=[   1   3 127] Mg
  flux global =  2.5535E-01, i_nu =   1, lambda =   3334.08 AA
"""


def _monitor_run(tmp_path, log=MONITOR_LOG):
    run = tmp_path / 'run'
    run.mkdir()
    (run / 'params.nml').write_text(
        "&rt_iteration_params maxiter=99 /\n&nlte_params convlim=1d-3 /\n")
    log_file = tmp_path / 'run.log'
    log_file.write_text(log)
    return run, log_file


def test_monitor_counts_iterations_not_flux_lines(tmp_path):
    """Each NLTE iteration reports its populations and its flux; only one of them is an
    iteration. Counting both doubled the progress and ran it past its own total."""
    from tumult.launch import monitor

    run, log_file = _monitor_run(tmp_path)
    done, total = monitor(run, log_file)
    assert (done, total) == (2, 99)         # two iterations after the LTE start, of maxiter


def test_monitor_reports_what_it_can_without_a_log(tmp_path):
    """A run launched with no log file, or read before it wrote one, must still answer."""
    from tumult.launch import monitor

    run, log_file = _monitor_run(tmp_path)
    assert monitor(run) == (0, 0)
    assert monitor(run, tmp_path / 'absent.log') == (0, 0)
    assert monitor(tmp_path / 'absent') == (0, 0)
    log_file.write_text('')                 # launched, nothing logged yet
    assert monitor(run, log_file) == (0, 99)


def test_progress_watch_reports_a_running_run(tmp_path):
    """The watcher turns the run's own files into a fraction while the work blocks."""
    from tumult.synthesis import progress

    run, log_file = _monitor_run(tmp_path, log=MONITOR_LOG.split('Iteration  2')[0])
    seen = []
    with progress.watch(run, lambda *args: seen.append(args), log_file=log_file,
                        interval=0.02) as watcher:
        assert watcher.poll() and watcher.fraction == pytest.approx(1 / 99)
        log_file.write_text(MONITOR_LOG)
        deadline = time.time() + 3
        while not seen and time.time() < deadline:
            time.sleep(0.02)
    assert seen[-1] == (2, 99, pytest.approx(2 / 99))


def test_progress_watch_survives_a_failing_callback(tmp_path):
    """A closed document must not take down the job whose progress it was showing."""
    from tumult.synthesis import progress

    run, log_file = _monitor_run(tmp_path, log=MONITOR_LOG.split('Iteration  1')[0])
    calls = []

    def report(*args):
        calls.append(args)
        raise RuntimeError('document gone')

    with progress.watch(run, report, log_file=log_file, interval=0.02) as watcher:
        log_file.write_text(MONITOR_LOG)
        deadline = time.time() + 3
        while not calls and time.time() < deadline:
            time.sleep(0.02)
    assert calls and watcher.done == 2


def test_snapshot_phases_split_the_settings():
    """The departures phase drops the spectrum; the spectrum phase drops the solve."""
    from tumult.synthesis import snapshot

    settings = {'atom_params': {'atom_file': 'mg.txt'},
                'linelist_params': {'linelist_file': ['master.txt']},
                'spectrum_params': {'R': 60000},
                'rt_iteration_params': {'maxiter': 3}}

    dep = snapshot.departure_settings(settings, '/snap', (30, 30, 128), 'MUST')
    assert 'linelist_params' not in dep and 'spectrum_params' not in dep
    assert dep['rt_output_params']['angle_quad'] == 'none'   # no formal solution
    # Explicit solver limits survive phase planning; 99 is only the default when absent.
    assert dep['rt_iteration_params']['maxiter'] == 3
    assert dep['atmos_params'] == {'atmos_file': '/snap', 'save_atmos': True,
                                   'atmos_format': 'MUST', 'nx': 30, 'ny': 30, 'nz': 128}
    assert settings['rt_iteration_params']['maxiter'] == 3   # the caller's dict is untouched

    field = snapshot.DepartureField(folder=Path('/x/data/snap/dep_ab'),
                                    name=Path('snap/dep_ab'), snapshot='/snap',
                                    resolution=(30, 30, 128))
    spec = snapshot.spectrum_settings(field, settings, (5160, 5190))
    assert spec['rt_iteration_params']['maxiter'] == 0       # populations frozen
    assert spec['atmos_params'] == {'atmos_format': 'save', 'atmos_file': '/x/data/snap/dep_ab',
                                    'read_departures': True}
    assert spec['linelist_params'] == settings['linelist_params']
    assert (spec['spectrum_params']['aa_blue'], spec['spectrum_params']['aa_red']) == (5160., 5190.)


def test_snapshot_spectrum_identity_separates_windows():
    """Two windows off one field are two runs; _dep_bucket's digest cannot tell them apart."""
    from tumult.synthesis import snapshot
    from tumult.synthesis.api import _dep_bucket

    field = snapshot.DepartureField(folder=Path('/x/data/snap/dep_ab'),
                                    name=Path('snap/dep_ab'), snapshot='/snap',
                                    resolution=(30, 30, 128))
    blue = snapshot.spectrum_settings(field, {}, (5160, 5190))
    red = snapshot.spectrum_settings(field, {}, (5520, 5540))
    assert snapshot._spectrum_name(blue) != snapshot._spectrum_name(red)
    assert _dep_bucket('snap', blue) == _dep_bucket('snap', red)   # why the digest differs


def test_snapshot_departures_are_reused_across_windows(tmp_path, monkeypatch):
    """One solve, two windows: the second window must not re-enter the departures phase."""
    import tumult as tm
    from tumult.synthesis import snapshot

    execution = tmp_path / 'execution'
    execution.mkdir()
    launched = execution / 'launched.txt'
    _stub_dispatch(execution,
                   'with open("launched.txt", "a") as f: f.write(sys.argv[1]+"\\n")\n'
                   '(out/"save").mkdir(exist_ok=True)\n'
                   '(out/"save"/"pop_Mg_meta.txt").write_text("")\n'
                   '(out/"abund.txt").write_text("Mg 7.58\\n")\n')
    monkeypatch.chdir(execution)
    # the run reader is not under test: a stub dispatch.x writes no readable output
    monkeypatch.setattr(tm, 'read', lambda path: SimpleNamespace(folder=path))

    field = snapshot.solve_departures('snap3d', '/snap', {'atom_params': {'atom_file': 'mg'}},
                                      resolution=(3, 3, 8), atmos_format='MUST')
    assert field.species == ('Mg',) and field.resolution == (3, 3, 8)
    assert field.folder == execution / 'data' / field.name

    snapshot.synthesise(field, window=(5160, 5190))
    snapshot.synthesise(field, window=(5520, 5540))
    snapshot.synthesise(field, window=(5160, 5190))           # a repeat, served from disk

    names = [line for line in launched.read_text().split() if line]
    assert sum('dep_' in n for n in names) == 1               # one solve
    assert len({n for n in names if 'spec_' in n}) == 2       # two distinct spectra
    assert sum('spec_' in n for n in names) == 2              # the repeat did not re-run


# --- durable workflow handles ------------------------------------------------------------
def _stub_request(root):
    """A resolvable request whose execution is stubbed; nothing here runs dispatch."""
    from tumult.synthesis.request import SynthRequest
    return SynthRequest(prefs={'abund_file': None, 'dispatch_cwd': str(root)},
                        grid={'folder': ''},
                        run_vals={'abundances': '', 'relative_abundances': False,
                                  'mode': 'LTE'},
                        target=(5750.0, 4.4, 0.0, 1.0),
                        masks=({'fit': (5000., 5001.), 'synth': (5000., 5001.)},))


def _stub_result(resolved):
    import numpy as np
    from tumult.spectrum import Spectrum
    from tumult.synthesis.request import SynthesisResult
    xx = np.linspace(5000., 5001., 5)
    return SynthesisResult((Spectrum(xx, np.ones_like(xx)),), (), resolved.request, {})


def test_workflow_cancelled_while_queued_settles_without_a_worker(tmp_path, monkeypatch):
    """Cancelling a queued workflow cancels its future, so no pool thread is ever spent
    on it, and its record still ends as cancelled; the live handle table lets it go."""
    from tumult.synthesis import submission, workflows

    ran = []
    monkeypatch.setattr(workflows, 'execute_request',
                        lambda resolved, cancel_token=None: ran.append(1) or _stub_result(resolved))
    gate = threading.Event()
    pool = concurrent.futures.ThreadPoolExecutor(1)
    monkeypatch.setattr(submission, 'WORKFLOWS', pool)
    blocker = pool.submit(gate.wait, 5)
    try:
        handle = submission.submit(_stub_request(tmp_path))
        seen = []
        handle.add_done_callback(lambda h: seen.append(h.status))
        assert handle.status == 'queued' and handle.id in submission._live
        assert handle.cancel() is True
        assert seen == ['cancelled']        # the callback sees the settled record
        assert handle.done and handle.status == 'cancelled'
        handle.add_done_callback(lambda h: seen.append(h.status))   # late registration
        assert seen == ['cancelled'] * 2
        assert handle.id not in submission._live
        with pytest.raises(CancelledError):
            handle.result()
        assert handle.cancel() is False              # already settled
    finally:
        gate.set()
        blocker.result(3)
        pool.shutdown(wait=True)
    assert not ran
    reopened = submission.open_workflow(handle.folder)
    assert reopened.status == 'cancelled' and reopened.done


def test_workflow_cancelled_while_running_and_settled_result_is_reloadable(tmp_path, monkeypatch):
    from tumult.synthesis import submission, workflows

    started, release = threading.Event(), threading.Event()

    def execute(resolved, cancel_token=None):
        started.set()
        assert release.wait(5)
        cancel_token.check()
        return _stub_result(resolved)

    monkeypatch.setattr(workflows, 'execute_request', execute)
    handle = submission.submit(_stub_request(tmp_path))
    assert started.wait(5) and handle.status == 'running'
    assert handle.cancel() is True
    release.set()
    with pytest.raises(CancelledError):
        handle.result(timeout=5)
    assert handle.status == 'cancelled' and handle.id not in submission._live

    started.clear()
    release.set()
    handle = submission.submit(_stub_request(tmp_path))
    result = handle.result(timeout=5)
    assert handle.status == 'finished' and len(result.spectra) == 1
    assert submission.open_workflow(handle.folder).result().spectra[0].values[0][0] == 5000.


def test_workflow_submission_failure_is_recorded(tmp_path, monkeypatch):
    from tumult.synthesis import submission

    pool = concurrent.futures.ThreadPoolExecutor(1)
    pool.shutdown()
    monkeypatch.setattr(submission, 'WORKFLOWS', pool)
    with pytest.raises(RuntimeError):
        submission.submit(_stub_request(tmp_path))
    handles = submission.discover_workflows(tmp_path / 'synth')
    assert [h.status for h in handles] == ['failed']
    assert handles[0].id not in submission._live
    assert 'RuntimeError' in handles[0].to_dict()['error']


# --- termination causes of a cancellable native run -----------------------------------
def test_cancellable_launch_records_success_failure_and_cancellation(tmp_path, monkeypatch):
    """The token-managed launch branch keeps the original cause: a clean exit is
    finished, a nonzero status is failed with its exit code and log, and only an
    explicit cancellation is recorded as cancelled (the killed child's status is not)."""
    from tumult import launch
    from tumult.synthesis import jobrecord

    execution = tmp_path / 'execution'
    execution.mkdir()
    monkeypatch.chdir(execution)
    _stub_dispatch(execution, '(out/"abund.txt").write_text("done")\n')
    output, _ = launch('cache/ok', verbose=0, cancel_token=CancelToken())
    assert jobrecord.state_of(output) == 'finished'

    _stub_dispatch(execution, '(out/"abund.txt").write_text("partial")\nsys.exit(3)\n')
    token = CancelToken()
    with pytest.raises(RuntimeError, match='status 3'):
        launch('cache/crash', verbose=0, cancel_token=token)
    record = jobrecord.read(execution / 'data/cache/crash')
    assert jobrecord.state_of(execution / 'data/cache/crash') == 'failed'
    assert record['returncode'] == 3 and record['error'].startswith('RuntimeError')
    assert Path(record['log']).name == 'crash.log' and Path(record['log']).is_file()
    assert not token.cancelled                     # a crash is not a cancellation
    assert not (execution / 'data/cache/crash/abund.txt').exists()

    _stub_dispatch(execution, 'time.sleep(30)\n')
    token = CancelToken()
    threading.Timer(0.4, token.cancel).start()
    with pytest.raises(CancelledError):
        launch('cache/stopped', verbose=0, cancel_token=token)
    record = jobrecord.read(execution / 'data/cache/stopped')
    assert record['state'] == 'cancelled' and record['returncode'] < 0
    assert record['error'] is None


# --- bounded server shutdown ------------------------------------------------------------
def test_shutdown_is_bounded_and_reports_uninterruptible_work(tmp_path, monkeypatch):
    """Shutdown stops admission, cancels queued workflows, kills a cancellable dispatch
    child and returns within its bound, naming the native worker that ignored it."""
    from tumult import launch
    from tumult.synthesis import runtime, submission, workflows
    from tumult.synthesis import jobrecord

    pools = {name: concurrent.futures.ThreadPoolExecutor(2, thread_name_prefix=f'test-{name}')
             for name in ('WORKFLOWS', 'HOP', 'IO')}
    for name, pool in pools.items():
        monkeypatch.setattr(runtime, name, pool)
    monkeypatch.setattr(submission, 'WORKFLOWS', pools['WORKFLOWS'])
    monkeypatch.setattr(runtime, '_default', ExecutionService(2))

    # 1. an in-process native call that cannot be interrupted
    release = threading.Event()
    native = pools['HOP'].submit(release.wait, 30)

    # 2. a running workflow honouring its token, and one queued behind it
    started = threading.Event()

    def execute(resolved, cancel_token=None):
        started.set()
        while not cancel_token.cancelled:
            time.sleep(0.01)
        cancel_token.check()

    monkeypatch.setattr(workflows, 'execute_request', execute)
    pools['WORKFLOWS'].submit(lambda: None).result(5)
    gate = threading.Event()
    blocker = pools['WORKFLOWS'].submit(gate.wait, 30)      # occupies the second thread
    running = submission.submit(_stub_request(tmp_path))
    assert started.wait(5)
    queued = submission.submit(_stub_request(tmp_path))
    assert queued.status == 'queued'

    # 3. a cancellable dispatch child on the leaf service
    execution = tmp_path / 'execution'
    execution.mkdir()
    ready = tmp_path / 'ready'
    _stub_dispatch(execution, f'Path({str(ready)!r}).touch()\ntime.sleep(30)\n')

    def run(token, slots):
        return launch('cache/child', verbose=0, cwd=execution, cancel_token=token)

    child = runtime.get_runtime().submit('child', 'sig', 'owner', run, token=CancelToken())
    deadline = time.monotonic() + 5
    while not ready.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    assert ready.exists()

    t0 = time.monotonic()
    report = runtime.shutdown(timeout=2.0)
    elapsed = time.monotonic() - t0
    assert elapsed < 4.0, elapsed
    assert report['stopped'] is False
    # the two threads that ignored everything: the HOP call and the pool blocker
    assert sorted(name.split(':')[0] for name in report['running']) == ['hop', 'workflows']
    assert report['cancelled_workflows'] == 2
    assert queued.status == 'cancelled' and running.status == 'cancelled'
    with pytest.raises(CancelledError):
        child.result(1)
    assert jobrecord.state_of(execution / 'data/cache/child') == 'cancelled'
    with pytest.raises(RuntimeError):
        submission.submit(_stub_request(tmp_path))        # admission is closed
    with pytest.raises(RuntimeError):
        runtime.get_runtime().submit('late', 'sig', 'owner', run)
    assert not native.done()                               # never claimed stopped
    release.set()
    gate.set()
    native.result(5)
    blocker.result(5)


def test_discovery_is_scoped_to_the_workflow_folder_and_keeps_the_label(tmp_path, monkeypatch):
    """``discover_workflows`` lists only ``<root>/_workflows/*``, newest first, with the
    submitter's label and meta; a record left running by a vanished process reports
    ``interrupted`` and cannot be cancelled (no PID is ever signalled)."""
    import json
    from tumult.synthesis import submission, workflows

    monkeypatch.setattr(workflows, 'execute_request',
                        lambda resolved, cancel_token=None: _stub_result(resolved))
    handle = submission.submit(_stub_request(tmp_path), label='t5750', meta={'k': [1]})
    handle.result(timeout=5)
    assert handle.label == 't5750' and handle.meta == {'k': [1]}

    stale = tmp_path / 'synth' / '_workflows' / '20260901T000000-stale'
    stale.mkdir(parents=True)
    (stale / submission.RECORD).write_text(json.dumps(
        {'version': 1, 'id': stale.name, 'state': 'running', 'phase': 'synthesis',
         'submitted': 1.0, 'request': {}, 'runs': [], 'result': None}))
    elsewhere = tmp_path / 'synth' / 'cache' / 'run' / '_workflows' / 'x'
    elsewhere.mkdir(parents=True)
    (elsewhere / submission.RECORD).write_text('{"state": "finished", "submitted": 9e9}')

    found = submission.discover_workflows(tmp_path / 'synth')
    assert [h.id for h in found] == [handle.id, stale.name]
    assert found[1].status == 'interrupted' and found[1].label == stale.name
    assert found[1].done and found[1].cancel() is False
    with pytest.raises(RuntimeError, match='interrupted'):
        found[1].result()
    assert [h.id for h in submission.discover_workflows(tmp_path / 'synth', limit=1)] == [handle.id]
