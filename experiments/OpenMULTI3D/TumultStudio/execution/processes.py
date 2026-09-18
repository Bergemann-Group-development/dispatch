"""The Processes panel, the job details modal and the recovery modal, rendered as HTML.

A panel row has two lines: what the job is (``Kind · name``) and how it runs (``8 threads``,
``queued``, ``stopping…``, ``failed · <error>``, ``Slurm <id>``), followed by the elapsed
time, which the page fills in from ``data-since`` once a job has run for 10 s. The rows
are this document's jobs (:class:`jobs.JobTable`), then the daemon's queued and running
workflows this document does not follow (other tabs, scripts) and the live dispatch runs
started outside any workflow, both from :class:`tumult.synthesis.daemon.LiveJobs`. A job
admitted less than ``SHOW_AFTER_S`` ago is not drawn, so work that settles at once (a
cached load) never flickers a row.

Server-side Python jobs (loads, line strengths, fits of broadening, normalisations, grid
transfers) run on one thread each and are not counted in the header, which shows the
daemon's threads in use plus those of external runs on this host over its budget.
"""
from __future__ import annotations

import html
import json
from pathlib import Path
import time

from .jobs import CANCELLING, FAILED, FINISHED, QUEUED, RUNNING, ExternalRun

SHOW_AFTER_S = 0.5
LOG_TAIL_LINES = 20

# the kind of a job, from its operation or the prefix of its internal name
KIND_LABELS = {
    "synthesis.synthesize": "Synthesis", "synthesis.fit": "Fit",
    "abundance.nlte_corrections": "NLTE corrections", "abundance.fit": "Abundance synthesis",
    "sparams.ionbalance": "Ionisation balance", "departures.solve": "Departure solve",
    "dispatch": "Dispatch run",
}
NAME_PREFIXES = {
    "abund-lines": "Line measurement", "rungs": "Abundance rungs",
    "linemask-auto": "Line mask", "import": "Import", "list": "Grid listing",
    "download": "Grid download",
}
NAME_LABELS = {
    "line-strengths": "Line strengths", "broadening-fit": "Broadening fit",
    "normalise-peaks": "Normalisation", "normalise-selection": "Normalisation",
}


def job_title(name, kind=None, path=None):
    """``Kind · name`` for a job's internal name and, for a durable workflow, its task
    kind. A durable workflow's name is the name of the object it delivers to
    (:mod:`naming`); a server-side job's name may carry a kind prefix
    (``abund-lines:<result>``, ``download MARCS``)."""
    name = str(name or "")
    if name in NAME_LABELS:
        return NAME_LABELS[name]
    if kind is not None:
        return f"{KIND_LABELS.get(kind, kind)} · {name}" if name else KIND_LABELS.get(kind, kind)
    rest, label = name, None
    for separator in (":", " "):
        prefix, found, tail = name.partition(separator)
        if found and prefix in NAME_PREFIXES:
            label, rest = NAME_PREFIXES[prefix], tail
            break
    if label is None and path is not None:
        label = "Load"
    return f"{label} · {rest}" if label and rest else (label or rest)


def _threads(n):
    return f"{n} thread{'' if n == 1 else 's'}"


def _elapsed(seconds):
    seconds = int(max(0.0, seconds))
    if seconds < 60:
        return f"{seconds} s"
    if seconds < 3600:
        return f"{seconds // 60} min"
    return f"{seconds // 3600} h {seconds % 3600 // 60:02d} min"


def _button(action, uid, icon, title, danger=False, confirm=False):
    return (f'<button class="object-btn{" danger" if danger else ""}" data-action="{action}" '
            f'data-uid="{html.escape(uid)}" data-confirm="{1 if confirm else 0}" '
            f'title="{html.escape(title)}"><i class="ti {icon}"></i></button>')


def _row(uid, title, line, state, since=None, buttons=(), error=None):
    since_attr = f' data-since="{since:.3f}"' if since else ""
    icon = ('<i class="ti ti-alert-circle"></i>' if state == FAILED
            else '<span class="spinner"></span>')
    tip = html.escape(error or line)
    return (f'<div class="object-item process{" failed" if state == FAILED else ""}" '
            f'data-pid="{html.escape(uid)}" data-state="{state}">{icon}'
            f'<div class="process-meta"><span class="object-name" title="{html.escape(title)}">'
            f'{html.escape(title)}</span>'
            f'<span class="process-state" title="{tip}">{html.escape(line)}'
            f'<span class="process-elapsed"{since_attr}></span></span></div>'
            f'{"".join(buttons)}</div>')


def _job_row(job, live):
    """A row for one job of this document."""
    summary = live.get(job.workflow_id) if live is not None and job.handle is not None else None
    state = job.display_state(summary)
    kind = (summary.get("kind") if summary else None) or job.operation
    if job.handle is not None and kind is None:
        try:
            kind = job.handle.kind.name
        except (OSError, KeyError, ValueError):
            kind = None
    title = job_title(job.name, kind, job.path)
    details = job.handle is not None or state == FAILED
    buttons = [_button("job-details", f"job-{job.pid}", "ti-info-circle", "Details")] if details else []
    if state == FAILED:
        buttons.append(_button("dismiss-job", str(job.pid), "ti-x", "Dismiss"))
        error = job.error or "failed"
        return _row(f"job-{job.pid}", title, f"failed · {error}", FAILED, buttons=buttons,
                    error=error)
    if state != CANCELLING:
        buttons.append(_button("cancel-pending", str(job.pid), "ti-x", "Cancel", danger=True,
                               confirm=job.heavy))
    since = None
    if state == QUEUED:
        line = "queued"
    elif state == CANCELLING:
        line = "stopping…"
    elif job.handle is not None:
        driver = (summary or {}).get("driver") or job.handle.detached
        started = (summary or {}).get("started") or job.started
        since = started
        if driver == "slurm":
            slurm_id = (summary or {}).get("slurm_id") or (
                (job.handle.to_dict().get("driver") or {}).get("slurm_id"))
            line = f"Slurm {slurm_id}" if slurm_id else "Slurm"
        else:
            threads = (summary or {}).get("threads")
            line = _threads(threads) if threads else "running"
    else:
        line, since = _threads(1), job.started
    return _row(f"job-{job.pid}", title, line, state, since, buttons)


def _live_row(summary):
    """A row for a daemon job this document does not follow: a workflow of another tab
    or a script, or an external dispatch run."""
    uid = summary["job"]
    external = summary["type"] == "external"
    title = job_title(summary.get("label") or uid, summary.get("kind"))
    state = summary.get("state") or RUNNING
    if state == "submitted" or state == QUEUED:
        line = "queued"
        state = QUEUED
    elif summary.get("slurm_id"):
        line, state = f"Slurm {summary['slurm_id']}", RUNNING
    elif summary.get("threads"):
        line, state = _threads(summary["threads"]), RUNNING
    else:
        line, state = "running", RUNNING
    buttons = [_button("job-details", uid, "ti-info-circle", "Details"),
               _button("run-cancel" if external else "workflow-cancel", uid, "ti-x",
                       "End this run" if external else "Cancel", danger=True, confirm=True)]
    since = summary.get("started") if state == RUNNING else None
    return _row(uid, title, line, state, since, buttons)


def processes_html(jobs, live=None, recovery=(), now=None):
    """The Processes panel body (see the module docstring), then a link to the recovery
    modal when records wait there."""
    now = time.time() if now is None else now
    rows, followed = [], set()
    for job in jobs:
        if job.workflow_id:
            followed.add(job.workflow_id)
        age_from = job.started if job.handle is not None and job.started else job.submitted
        if not job.terminal and now - age_from < SHOW_AFTER_S:
            continue
        rows.append(_job_row(job, live))
    if live is not None:
        for summary in live.jobs():
            if summary["job"] not in followed:
                rows.append(_live_row(summary))
    count = len(recovery)
    if count:
        rows.append('<button class="object-btn process-recover" data-action="workflow-discover" '
                    f'data-uid="-">{count} record{"" if count == 1 else "s"} to recover</button>')
    return "".join(rows)


def usage_html(live):
    """The header's thread count: the daemon's threads in use plus those of external runs
    on this host, over the daemon's budget."""
    if live is None or not live.budget:
        return ""
    return f"{live.used + live.external_threads()}/{live.budget} threads"


# --- the details modal ----------------------------------------------------------------------
def _section(title, body):
    return f'<div class="job-section"><div class="job-section-title">{title}</div>{body}</div>' if body else ""


def _pairs(items):
    rows = "".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>"
                   for k, v in items if v not in (None, "", ()))
    return f'<table class="job-pairs">{rows}</table>' if rows else ""


def _when(epoch):
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(epoch)) if epoch else None


def _tail(path, lines=LOG_TAIL_LINES):
    try:
        with open(path, "rb") as handle:
            handle.seek(0, 2)
            size = handle.tell()
            handle.seek(max(0, size - 8192))
            text = handle.read().decode(errors="replace")
    except OSError:
        return ""
    return "\n".join(text.splitlines()[-lines:])


def _logs(paths):
    return "".join(
        f'<details class="job-log"><summary title="{html.escape(str(path))}">'
        f'{html.escape(Path(path).name)}</summary><pre>{html.escape(_tail(path))}</pre></details>'
        for path in paths)


def _inputs(meta):
    if not meta:
        return ""
    params = meta.get("params") if isinstance(meta, dict) else None
    if params:
        return _pairs((row[0], row[1]) for row in params if len(row) >= 2)
    return f"<pre>{html.escape(json.dumps(meta, indent=1, default=str))}</pre>"


def _runs_table(runs):
    if not runs:
        return ""
    rows = "".join(
        f'<tr><td title="{html.escape(run["output"])}">{html.escape(Path(run["output"]).name)}</td>'
        f'<td>{html.escape(run["state"])}</td><td>{run["progress"][0]}/{run["progress"][1]}</td></tr>'
        for run in runs)
    return f'<table class="job-runs"><tr><th>run</th><th>state</th><th>progress</th></tr>{rows}</table>'


def _actions(buttons):
    return f'<div class="job-actions">{"".join(buttons)}</div>'


def _text_button(action, uid, label, danger=False):
    return (f'<button class="text-btn job-action{" danger" if danger else ""}" '
            f'data-action="{action}" data-uid="{html.escape(uid)}">{html.escape(label)}</button>')


def _copy_button(value, label):
    return (f'<button class="text-btn job-copy" data-copy="{html.escape(value)}">'
            f'{html.escape(label)}</button>')


def workflow_details_html(handle, summary=None, job=None, recoverable=None):
    """The details of a durable workflow: progress, threads and driver, timing, inputs,
    its runs, logs and the actions that apply to it."""
    record = handle.to_dict()
    progress = handle.progress
    state = record.get("state")
    driver = record.get("driver") or {}
    diagnostics = progress.get("diagnostics") or {}
    title = job_title(record.get("label") or handle.id, (summary or {}).get("kind") or handle.kind.name)
    now = time.time()
    started, ended = record.get("started"), record.get("ended")
    body = [f'<div class="job-title">{html.escape(title)}</div>',
            _section("Progress", _pairs([
                ("state", state), ("phase", record.get("phase")),
                ("runs", f"{diagnostics.get('completed', 0)}/{diagnostics['native_runs']} finished"
                 if diagnostics.get("native_runs") else len(record.get("runs") or ()) or None),
                ("error", record.get("error"))])),
            _section("Execution", _pairs([
                ("threads", (summary or {}).get("threads") if state == RUNNING else None),
                ("driver", driver.get("kind")), ("pid", driver.get("pid")),
                ("host", driver.get("host")), ("Slurm job", driver.get("slurm_id"))])),
            _section("Timing", _pairs([
                ("submitted", _when(record.get("submitted"))), ("started", _when(started)),
                ("ended", _when(ended)),
                ("elapsed", _elapsed((ended or now) - started) if started else None)])),
            _section("Inputs", _inputs(record.get("meta"))),
            _section("Runs", _runs_table(diagnostics.get("runs"))),
            _section("Logs", _logs(handle.logs))]
    buttons = []
    if state in (QUEUED, RUNNING):
        if job is not None:
            buttons.append(_text_button("cancel-pending", str(job.pid), "Cancel", danger=True))
        else:
            buttons.append(_text_button("workflow-attach", handle.id, "Attach"))
            buttons.append(_text_button("workflow-cancel", handle.id, "Cancel", danger=True))
    if recoverable is not None:
        if state == FINISHED:
            buttons.append(_text_button("workflow-import", handle.id, "Import"))
        buttons.append(_text_button("workflow-forget", handle.id, "Forget"))
    buttons += [_copy_button(handle.id, "Copy id"),
                _text_button("folder-open", handle.id, "Open folder")]
    return "".join(body) + _actions(buttons)


def run_details_html(entry):
    """The details of an external dispatch run (:class:`jobs.ExternalRun`)."""
    body = [f'<div class="job-title">{html.escape(job_title(entry.name, "dispatch"))}</div>',
            _section("Execution", _pairs([
                ("state", entry.state), ("threads", entry.threads), ("pid", entry.pid),
                ("host", entry.host), ("Slurm job", entry.slurm_id), ("input", entry.input),
                ("folder", entry.folder)])),
            _section("Timing", _pairs([("started", _when(entry.started))])),
            _section("Logs", _logs(entry.logs))]
    buttons = []
    if entry.live:
        buttons.append(_text_button("run-cancel", entry.id, "End this run", danger=True))
    if entry.departures:
        buttons.append(_text_button("run-departures", entry.id, "Use as saved departures"))
    buttons += [_copy_button(entry.folder, "Copy folder"),
                _text_button("folder-open", entry.id, "Open folder")]
    return "".join(body) + _actions(buttons)


def job_details_html(job):
    """The details of a failed Python job of this document: its error and logs."""
    body = [f'<div class="job-title">{html.escape(job_title(job.name, path=job.path))}</div>',
            _section("Outcome", _pairs([("state", job.state), ("error", job.error),
                                        ("file", job.path)])),
            _section("Timing", _pairs([("submitted", _when(job.submitted)),
                                       ("ended", _when(job.ended))])),
            _section("Logs", _logs(job.logs))]
    return "".join(body) + _actions([_text_button("dismiss-job", str(job.pid), "Dismiss")])


# --- the recovery modal -----------------------------------------------------------------------
def recovery_html(recovery):
    """The settled workflow records this document has not imported, newest first: a
    finished one can be imported, any of them opened in details or forgotten."""
    entries = list(recovery)
    if not entries:
        return '<div class="file-empty">No workflow records to recover.</div>'
    rows = []
    for entry in entries:
        wid = entry.id
        if entry.state == FINISHED:
            icon, note = "ti-package", "finished · result on disk"
        elif entry.state == "interrupted":
            icon, note = "ti-plug-x", "interrupted: the process that ran it is gone; nothing resumes"
        elif entry.live:
            icon, note = "ti-plug", f"{entry.state} · not attached"
        else:
            icon = "ti-alert-circle" if entry.state == FAILED else "ti-ban"
            note = entry.state + (f": {entry.error}" if entry.error else "")
        if entry.referenced:
            note += " · from the loaded save"
        buttons = [_button("job-details", wid, "ti-info-circle", "Details")]
        if entry.live:
            buttons.append(_button("workflow-attach", wid, "ti-plug-connected",
                                   "Attach: follow it here and allow cancelling"))
        elif entry.state == FINISHED:
            buttons.append(_button("workflow-import", wid, "ti-download",
                                   "Import its result into the workspace"))
        buttons.append(_button("workflow-forget", wid, "ti-x",
                               "Forget this record here (it stays on disk)"))
        logs = "".join(f'<span class="process-log" title="{html.escape(p)}">'
                       f'{html.escape(Path(p).name)}</span>' for p in entry.logs)
        rows.append(
            f'<div class="object-item process recoverable" data-wid="{html.escape(wid)}" '
            f'data-state="{html.escape(entry.state)}"><i class="ti {icon}"></i>'
            f'<div class="process-meta"><span class="object-name" title="{html.escape(wid)}">'
            f'{html.escape(job_title(entry.name, entry.kind))}</span>'
            f'<span class="process-state" title="{html.escape(note)}">{html.escape(note)}</span>'
            f'{logs}</div>{"".join(buttons)}</div>')
    return "".join(rows)
