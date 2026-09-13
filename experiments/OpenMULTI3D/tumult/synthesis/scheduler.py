"""Thin Slurm client for detached workflows: submit the driver, poll it, cancel it.

Every call shells out to the site's ``sbatch``/``squeue``/``scancel``; nothing here
parses accounting output beyond the job id and the one-word state. A missing scheduler
raises ``FileNotFoundError`` from the first call, which the submitter reports.
"""
from __future__ import annotations

import re
import subprocess

_SUBMITTED = re.compile(r'Submitted batch job (\d+)')
# squeue reports these while the job still holds or awaits an allocation
_PENDING = {'PENDING', 'CONFIGURING', 'REQUEUED', 'REQUEUE_HOLD', 'SUSPENDED'}
_ACTIVE = {'RUNNING', 'COMPLETING', 'STAGE_OUT'}


def sbatch(script, args, *, env=None):
    """Submit ``script`` with positional ``args``; return the job id as a string."""
    result = subprocess.run(['sbatch', str(script), *map(str, args)], env=env,
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'sbatch failed: {result.stderr.strip() or result.stdout.strip()}')
    match = _SUBMITTED.search(result.stdout)
    if match is None:
        raise RuntimeError(f'sbatch gave no job id: {result.stdout.strip()}')
    return match.group(1)


def state(job_id):
    """``queued``, ``running`` or ``None`` once the scheduler no longer lists the job.
    ``squeue`` forgets a job shortly after it ends, so ``None`` says only that the
    driver is not going to write anything more; the record names the outcome."""
    result = subprocess.run(['squeue', '-h', '-j', str(job_id), '-o', '%T'],
                            capture_output=True, text=True)
    word = result.stdout.strip().split('\n')[0].strip() if result.stdout.strip() else ''
    if result.returncode or not word:
        return None
    if word in _ACTIVE:
        return 'running'
    if word in _PENDING:
        return 'queued'
    return None


def scancel(job_id):
    """Ask the scheduler to end the job; ``True`` when it accepted the request."""
    result = subprocess.run(['scancel', str(job_id)], capture_output=True, text=True)
    return result.returncode == 0
