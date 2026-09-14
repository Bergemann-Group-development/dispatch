"""Run one recorded workflow to completion on a Slurm node.

    python -m tumult.synthesis.driver <workflow folder>

This is the body of a ``detach='slurm'`` workflow: the process the batch script runs on
its allocated node. It reads the task from ``workflow.json``, records its pid and host
beside the driver block the submitter wrote, executes the task the way the workflow
daemon does and settles the record. SIGTERM, SIGINT and SIGHUP cancel its token, so
``scancel`` ends the run through the same path as a Processes-panel cancel: the dispatch
child is killed, its partial cache entry removed, and the record reads ``cancelled``.

Exit status: 0 finished, 1 failed, 2 cancelled, 3 no usable record.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
import signal
import socket
import sys
import time

from .cancel import CancelledError, CancelToken
from . import submission

log = logging.getLogger(__name__)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 3
    folder = Path(argv[0]).expanduser().resolve()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    record = submission._read(folder)
    if record is None:
        log.error('no %s in %s', submission.RECORD, folder)
        return 3
    if record.get('state') not in ('queued', 'running'):
        log.error('workflow %s is already %s', record.get('id'), record.get('state'))
        return 3
    try:
        kind, args = submission.decode_task(record)
    except Exception as exc:
        submission.update_record(folder, record, state='failed', phase='failed',
                                 ended=time.time(), error=f'unreadable request: {exc}')
        log.exception('workflow %s has an unreadable request', record.get('id'))
        return 3

    token = CancelToken()

    def on_signal(signum, _frame):
        log.info('workflow %s: signal %s, cancelling', record.get('id'),
                 signal.Signals(signum).name)
        token.cancel()

    for name in ('SIGTERM', 'SIGINT', 'SIGHUP'):
        if hasattr(signal, name):
            signal.signal(getattr(signal, name), on_signal)

    def claim(current):
        driver = current.get('driver') or {}
        driver.update(kind='slurm', pid=os.getpid(), host=socket.gethostname(),
                      started=time.time(),
                      slurm_id=driver.get('slurm_id') or os.environ.get('SLURM_JOB_ID'))
        current['driver'] = driver

    submission.mutate_record(folder, claim)
    log.info('workflow %s: driver pid %d on %s', record.get('id'), os.getpid(),
             socket.gethostname())
    try:
        submission.execute_recorded(folder, kind, args, token, record)
    except CancelledError:
        return 2
    except BaseException:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
