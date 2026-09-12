"""Document-owned result jobs with isolated inputs and coalesced progress."""
from copy import deepcopy
import logging
from dataclasses import dataclass
from threading import Lock

from tumult.synthesis.cancel import CancelToken, CancelledError


@dataclass
class _Outcome:
    result: object = None
    error: Exception | None = None


class ResultJobs:
    """One current submission per result; workers never mutate the displayed object."""

    def __init__(self, ctx, find=None, render=None):
        # ``find``/``render`` are only consulted for result-bound jobs; a job that
        # produces a new object (``result=None``) needs neither.
        self.ctx, self.find, self.render = ctx, find, render
        self.current = {}

    def invalidate(self, uid):
        self.current.pop(uid, None)

    def submit(self, name, work, done, *, result=None, on_error=None, detached=False):
        """Submit one job. ``detached`` keeps it running when the browser tab closes,
        for work too long to babysit; its result still only reaches a live document."""
        key = result.uid if result is not None else name
        generation = object()
        self.current[key] = generation
        working = deepcopy(result)
        lock = Lock()
        pending = []
        finished = False

        def valid():
            return (not getattr(self.ctx, "closed", False)
                    and self.current.get(key) is generation
                    and (result is None or self.find(key) is result))

        def apply(snapshot):
            result.__dict__.update(snapshot.__dict__)
            self.render(result)

        def flush():
            with lock:
                snapshot = pending.pop() if pending else None
            if snapshot is not None and valid() and not finished:
                apply(snapshot)

        def progress(value):
            # Coalesce before copying; published arrays/spectra are independent of the worker.
            with lock:
                if pending:
                    return
                pending.append(deepcopy(value))
            self.ctx.schedule(flush)

        def stopped(status):
            nonlocal finished
            if finished or not valid():
                return
            finished = True
            self.current.pop(key, None)
            if result is not None:
                result.termination = status
                if hasattr(result, "converged"):
                    result.converged = False
                self.render(result)

        def run(token):
            token = token or CancelToken()
            detach = token.add_callback(lambda: self.ctx.schedule(lambda: stopped("cancelled")))
            try:
                token.check()
                return _Outcome(work(working, token, progress))
            except CancelledError:
                raise
            except Exception as exc:
                # the toast shows only the message; the traceback goes to the server log
                logging.exception("%s failed", name)
                return _Outcome(error=exc)
            finally:
                detach()

        def complete(outcome):
            nonlocal finished
            if not valid() or finished:
                return
            if outcome.error is not None:
                stopped("failed")
                if on_error is not None:
                    on_error(outcome.error)      # the submitter has its own error surface
                else:
                    self.ctx.notify(f"{name} failed: {outcome.error}")
                return
            finished = True
            self.current.pop(key, None)
            if result is not None:
                apply(outcome.result)
                done(result)
            else:
                done(outcome.result)

        future = self.ctx.submit(name, run, complete, detached=detached)
        if future is not None:
            future.add_done_callback(
                lambda f: self.ctx.schedule(lambda: stopped("cancelled")) if f.cancelled() else None)
        return future
