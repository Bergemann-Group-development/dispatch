"""One browser scenario over a live TumultStudio2 server.

A Bokeh server runs the real directory app in a subprocess with a temporary config and
workflow root and a stubbed synthesis execution (no dispatch, no HOP); headless Chromium
drives the page through its own buttons, rows and modals, so every step exercises the
template's delegated events and the document callbacks behind them. Covered in order:
submitting a synthesis from the modal, cancelling it from its Processes row, a failed
workflow's retained details and their dismissal, saving a slot, replacing the workspace
from another slot, reconnecting as a new session that finds the record the previous one
left running, and importing its finished result once.

The stubbed execution follows a plan the test writes: each call waits for a release
file, fails, or finishes at once. Nothing here validates physics; the spectra are flat.

Requirements: ``node`` with the ``playwright`` npm package and its Chromium (the ``npx``
cache or ``TUMULT_NODE_PATH``; see README_INSTALL.md), else the test skips.
"""
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

SERVER = r"""
import itertools, json, os, sys, time
from pathlib import Path
ROOT, PORT, CTRL = Path(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
sys.path[:0] = [str(ROOT), str(ROOT / "TumultStudio2")]
import numpy as np
from tumult import Spectrum
from tumult.synthesis import workflows
from tumult.synthesis.request import SynthesisResult
calls = itertools.count()


def execute(resolved, cancel_token=None):
    step = json.loads((CTRL / "plan.json").read_text())[next(calls)]
    deadline = time.monotonic() + 60
    while step.get("wait") and not (CTRL / step["wait"]).exists():
        cancel_token.check()
        assert time.monotonic() < deadline, "release file never appeared"
        time.sleep(0.05)
    cancel_token.check()
    if step.get("fail"):
        raise RuntimeError(step["fail"])
    xx = np.linspace(5000., 5001., 5)
    return SynthesisResult((Spectrum(xx, np.ones_like(xx)),), (), resolved.request, {})


workflows.execute_request = execute
from bokeh.application import Application
from bokeh.application.handlers.directory import DirectoryHandler
from bokeh.server.server import Server
app = Application(DirectoryHandler(filename=str(ROOT / "TumultStudio2")))
server = Server({"/TumultStudio2": app}, port=PORT, allow_websocket_origin=[f"localhost:{PORT}"])
server.start()
server.io_loop.start()
"""

DRIVER = r"""
const { chromium } = require('playwright');
const readline = require('readline');
(async () => {
  const browser = await chromium.launch();
  let page = null;
  const errors = [];
  const open = async () => {
    page = await browser.newPage();
    page.on('pageerror', (e) => errors.push(String(e)));
    page.on('dialog', (d) => d.accept());
  };
  const rl = readline.createInterface({ input: process.stdin });
  for await (const line of rl) {
    const step = JSON.parse(line);
    let result = null;
    try {
      if (step.op === 'goto') {
        if (page) await page.close();
        await open();
        await page.goto(step.url, { waitUntil: 'networkidle' });
        await page.waitForFunction(() => window.Bokeh && Bokeh.documents.length > 0
          && Bokeh.documents[0].get_model_by_name('object-actions'), null, { timeout: 30000 });
      } else if (step.op === 'click') {
        await page.click(step.selector, { timeout: 10000 });
      } else if (step.op === 'fill') {
        await page.fill(step.selector, step.value, { timeout: 10000 });
      } else if (step.op === 'dispatch') {
        await page.dispatchEvent(step.selector, step.event, undefined, { timeout: 10000 });
      } else if (step.op === 'wait') {
        await page.waitForFunction(step.js, null, { timeout: step.timeout || 15000 });
      } else if (step.op === 'eval') {
        result = await page.evaluate(step.js);
      }
      console.log(JSON.stringify({ ok: true, result, errors: errors.splice(0) }));
      if (step.op === 'quit') break;
    } catch (e) {
      console.log(JSON.stringify({ ok: false, error: String(e).slice(0, 800),
                                   errors: errors.splice(0) }));
    }
  }
  await browser.close();
})().catch((e) => { console.error(e); process.exit(1); });
"""


def _node_playwright():
    explicit = os.environ.get("TUMULT_NODE_PATH")
    candidates = [Path(explicit)] if explicit else []
    candidates += sorted(Path.home().glob(".npm/_npx/*/node_modules"))
    for base in candidates:
        if (base / "playwright" / "package.json").exists():
            return base
    return None


def _free_port():
    with socket.socket() as sock:
        sock.bind(("localhost", 0))
        return sock.getsockname()[1]


class _Browser:
    """One headless Chromium, driven one JSON step at a time over the driver's stdin."""

    def __init__(self, tmp_path, node, base):
        script = tmp_path / "driver.js"
        script.write_text(DRIVER)
        self.proc = subprocess.Popen([node, str(script)], stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True, env={**os.environ, "NODE_PATH": str(base)})
        self.page_errors = []

    def step(self, op, **kw):
        self.proc.stdin.write(json.dumps({"op": op, **kw}) + "\n")
        self.proc.stdin.flush()
        line = self.proc.stdout.readline()
        assert line, "the browser driver exited"
        reply = json.loads(line)
        self.page_errors += reply.get("errors", [])
        assert reply["ok"], reply.get("error")
        return reply["result"]

    def click(self, selector):
        self.step("click", selector=selector)

    def fill(self, selector, value):
        self.step("fill", selector=selector, value=value)

    def wait(self, js, timeout=15000):
        self.step("wait", js=js, timeout=timeout)

    def eval(self, js):
        return self.step("eval", js=js)

    def text(self, selector):
        return self.eval(f"document.querySelector({selector!r}).textContent")

    def count(self, selector):
        return self.eval(f"document.querySelectorAll({selector!r}).length")

    def close(self):
        try:
            self.step("quit")
        finally:
            try:
                self.proc.communicate(timeout=15)
            except subprocess.TimeoutExpired:
                self.proc.kill()


_GRID = {"id": "g1", "name": "Test grid", "folder": "", "pattern": "", "format": "Marcs"}


def _seed_config(tmp_path, monkeypatch):
    """The config the server starts from: a grid that resolves without files, the
    workflow root under ``tmp_path`` and a synthesis form ready to submit."""
    import store
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "config.json")
    store.set_section("synthesis", {
        "dispatch_cwd": str(tmp_path), "synth_output_dir": "synth", "linelist": "",
        "abund_file": "", "absmet": "", "atmos_grids": [dict(_GRID)],
        "workflow_runner": "in-process",     # the stubbed execution lives in the server
        "atmos_grid_id": "g1", "ref_grid_id": "g1"})
    store.set_section("synth_run", {"teff": 5750.0, "logg": 4.4, "feh": 0.0, "vmic": 1.0,
                                    "mode": "LTE", "abundances": "", "atoms": [],
                                    "relative_abundances": True})
    return store.CONFIG_PATH


@pytest.fixture
def server(tmp_path, monkeypatch):
    """The directory app on a free port in a subprocess, stubbed execution included."""
    config = _seed_config(tmp_path, monkeypatch)
    ctrl = tmp_path / "ctrl"
    ctrl.mkdir()
    (ctrl / "plan.json").write_text("[]")
    port = _free_port()
    launcher = tmp_path / "server.py"
    launcher.write_text(SERVER)
    log = (tmp_path / "server.log").open("w")
    proc = subprocess.Popen([sys.executable, str(launcher), str(ROOT), str(port), str(ctrl)],
                            stdout=log, stderr=subprocess.STDOUT,
                            env={**os.environ, "TUMULT_CONFIG": str(config),
                                 "MPLCONFIGDIR": os.environ.get("MPLCONFIGDIR", "/tmp/tumult-mpl")})
    url = f"http://localhost:{port}/TumultStudio2"
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            pytest.fail("server exited: " + (tmp_path / "server.log").read_text()[-2000:])
        try:
            if urllib.request.urlopen(url, timeout=2).status == 200:
                break
        except OSError:
            time.sleep(0.3)
    else:
        pytest.fail("server never answered")
    try:
        yield url, ctrl, config
    finally:
        proc.terminate()
        proc.wait(15)
        log.close()


@pytest.fixture
def browser(tmp_path):
    node, base = shutil.which("node"), _node_playwright()
    if node is None or base is None:
        pytest.skip("node with a playwright package (npx cache or TUMULT_NODE_PATH) is "
                    "required for the browser scenario")
    driver = _Browser(tmp_path, node, base)
    try:
        yield driver
    finally:
        driver.close()


ROWS = "document.querySelectorAll('#processes .process:not(.recoverable)').length"
OBJECTS = "document.querySelectorAll('#objects .object-item').length"


def _submit(browser, xmin, xmax):
    """Open the modal from the landing page, set the window and press Synthesise. Each
    call needs a window the fields do not already show: the form only re-renders when its
    markup changes, and an unchanged input fires no change event to carry the value."""
    browser.click("#empty-synth")
    browser.fill('num-input[data-form="synthwindow"][data-field="xmin"] input', str(xmin))
    browser.fill('num-input[data-form="synthwindow"][data-field="xmax"] input', str(xmax))
    browser.click("#synth-modal .synth-submit")
    browser.wait("!document.getElementById('synth-modal').open")


def _records(tmp_path):
    """``{workflow id: recorded state}`` read from the records the server wrote (a
    handle opened here would call the server's live workflow interrupted)."""
    return {path.parent.name: json.loads(path.read_text())["state"]
            for path in (tmp_path / "synth" / "_workflows").glob("*/workflow.json")}



def _load_workspace(browser, panel, sid):
    """Select a saved workspace and load it; discard when the unsaved-changes prompt opens."""
    browser.click(f'{panel} .session-item[data-sid="{sid}"]')
    browser.click(f'{panel} .session-load')
    browser.eval("(function () { var d = document.getElementById('workspace-load-modal');"
                 " if (d.open) d.querySelector('[data-workspace-load=\"load-discard\"]').click();"
                 " return d.open; })()")


def test_browser_scenario_submit_cancel_fail_save_replace_reconnect_import(
        tmp_path, server, browser):
    url, ctrl, config = server
    (ctrl / "plan.json").write_text(json.dumps([
        {"wait": "release0"}, {"fail": "stub exploded"}, {"wait": "release2"}]))
    browser.step("goto", url=url)

    # submit and cancel: the row spins with its cancel button; the (confirmed) cancel
    # reaches the stub's token and the row leaves; the record on disk says cancelled
    _submit(browser, 5000, 5001)
    browser.wait(f"{ROWS} === 1 && document.querySelector('#processes .process .spinner')")
    assert "t5750g4.4m+0.0v1.0" in browser.text("#processes")
    browser.click('#processes [data-action="cancel-pending"]')
    browser.wait(f"{ROWS} === 0", timeout=20000)
    assert list(_records(tmp_path).values()) == ["cancelled"]
    assert browser.count("#objects .object-item") == 0

    # a failed workflow is retained with its error until dismissed; the footer points at
    # the panel
    _submit(browser, 5002, 5003)
    browser.wait("document.querySelector('#processes .process.failed')")
    panel = browser.text("#processes")
    assert "stub exploded" in panel and "t5750g4.4m+0.0v1.0" in panel
    assert "Processes" in browser.text("#footer-status")
    browser.click('#processes [data-action="dismiss-job"]')
    browser.wait(f"{ROWS} === 0")
    assert sorted(_records(tmp_path).values()) == ["cancelled", "failed"]

    # a third synthesis keeps running; save a slot that references it (the document
    # references every workflow it submitted, none imported yet)
    _submit(browser, 5004, 5005)
    browser.wait(f"{ROWS} === 1 && document.querySelector('#processes .process .spinner')")
    browser.click("#open-session")
    browser.fill("#session-name", "running")
    browser.click("#session-save")
    browser.wait("document.getElementById('footer-status').textContent.includes('Saved')")
    browser.click("#close-session")
    slot, = json.loads(config.read_text())["sessions"]
    running_id, = [wid for wid, state in _records(tmp_path).items() if state == "running"]
    assert slot["name"] == "running"
    assert {w["id"] for w in slot["workflows"]} == set(_records(tmp_path)) and running_id
    assert not any(w["imported"] for w in slot["workflows"])

    # open a spectrum through the picker, save that workspace too, then replace it with
    # the "running" slot: the object goes, the followed workflow's row stays
    spectrum = tmp_path / "obs.txt"
    spectrum.write_text("".join(f"{5000 + 0.01 * i:.3f} 1.000\n" for i in range(300)))
    browser.click("#open-file")
    browser.click('[data-import="spectrum"]')
    browser.fill('#file-modal input[type="text"]', str(spectrum))
    browser.step("dispatch", selector='#file-modal input[type="text"]', event="input")
    browser.wait(f"document.querySelector('#file-modal .file-item[data-path=\"{spectrum}\"]')")
    browser.click(f'#file-modal .file-item[data-path="{spectrum}"]')   # selects the row
    browser.click("#file-modal .picker-open")                             # opens it
    browser.wait(f"{OBJECTS} === 1")
    browser.click("#close-file")
    browser.click("#open-session")
    browser.fill("#session-name", "obs")
    browser.click("#session-save")
    browser.wait("document.querySelectorAll('.session-item').length >= 2")
    _load_workspace(browser, '#session-modal', 1)
    browser.wait(f"{OBJECTS} === 0")
    assert "Loaded" in browser.text("#footer-status")
    assert browser.count("#processes .process:not(.recoverable)") == 1

    # reconnect: a new session finds the record still running; once it finishes on the
    # server, loading the "obs" slot and pressing Recover offers the result for import
    browser.step("goto", url=url)
    browser.wait("document.querySelector('#processes .recoverable[data-state=\"running\"]')")
    assert browser.count('#processes [data-action="workflow-attach"]') == 1
    (ctrl / "release2").write_text("")
    deadline = time.monotonic() + 20
    while _records(tmp_path)[running_id] != "finished" and time.monotonic() < deadline:
        time.sleep(0.2)
    assert _records(tmp_path)[running_id] == "finished"
    _load_workspace(browser, '#canvas-empty', 2)
    browser.wait(f"{OBJECTS} === 1")
    browser.click('[data-action="workflow-discover"]')
    browser.wait("document.querySelector('#processes .recoverable[data-state=\"finished\"]')")
    browser.click('#processes [data-action="workflow-import"]')
    browser.wait(f"{OBJECTS} === 2")
    names = browser.eval("Array.from(document.querySelectorAll('#objects .object-name'))"
                         ".map(e => e.textContent)")
    assert names == ["obs.txt", "t5750g4.4m+0.0v1.0"]
    # once imported, the record is neither offered again nor importable: a second scan
    # (its footer report is a new message) lists no import button
    before = browser.eval("document.getElementById('footer-status').innerHTML")
    browser.click('[data-action="workflow-discover"]')
    browser.wait("document.getElementById('footer-status').innerHTML !== "
                 + json.dumps(before) + " && document.getElementById('footer-status')"
                 ".textContent.includes('recoverable')")
    assert browser.count('#processes [data-action="workflow-import"]') == 0
    assert browser.count("#objects .object-item") == 2
    assert not browser.page_errors
