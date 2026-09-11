"""Execution invariants: shared budgets, cancellation, cache writers and processes."""
import concurrent.futures
import os
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
