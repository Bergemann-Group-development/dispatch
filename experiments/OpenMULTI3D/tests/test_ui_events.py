"""The shared analysis-tool event wiring and the shared EW refresh.

Two halves. The headless half drives the real document's action channels and the
``filter-state`` event with a stubbed ``measure_ew``: what a refresh forwards, what it
assigns and what it leaves alone. The browser half renders the tools' actual markup into a
headless Chromium, lets the page's delegated listeners turn clicks into channel payloads
(``Bokeh.documents[0].get_model_by_name(...).setv`` is stubbed to record them), and feeds
every recorded payload back into the same document's named channels, so the DOM-side
routing and the server-side callbacks are exercised as one path. The browser is the one
node's ``playwright`` package drives; without it that half is skipped, not faked.
"""
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from tumult import Spectrum
from tumult import abundance as ab

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "TumultStudio2" / "templates" / "index.html"


@pytest.fixture(autouse=True)
def _temporary_config(tmp_path, monkeypatch):
    import store
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "tumultstudio.json")


def _build_document():
    from bokeh.document import Document
    from bokeh.io.state import curstate

    document = Document()
    curstate().document = document
    import app
    return document, app.build()


def _line(lo, hi, **fields):
    from tumult.abundance import LineResult
    from tumult.linemask import LineMask
    lam0 = 0.5 * (lo + hi) if np.isfinite(lo) else 5005.0
    mask = LineMask(lam0=lam0)
    if np.isfinite(lo):
        mask.set_fit(lo, hi, 0.5)              # else: a no-window line, nothing to refresh
    values = dict(lam0=lam0, ion="Fe I", e_low=1.0, loggf=-1.0, ew=50.0, ew_err=2.0,
                  mask=mask)
    return LineResult(**(values | fields))


def _stub_measure(monkeypatch, fail_at=()):
    """Replace ``measure_ew`` with a recorder: the width is the window's lower edge, the
    error the SNR (so both are recognisable), and ``fail_at`` windows raise."""
    calls = []

    def measure(spec, window=None, snr=None, cont_pct=0.0, profile="trapz", x0=None):
        calls.append(dict(spec=spec, window=tuple(window), snr=snr, cont_pct=cont_pct,
                          profile=profile, x0=x0))
        if window[0] in fail_at:
            raise ValueError("no profile")
        return SimpleNamespace(maa=float(window[0]), width_err_maa=float(snr or -1.0))
    monkeypatch.setattr(ab, "measure_ew", measure)
    return calls


def _observation(ctx, tmp_path):
    path = tmp_path / "obs.txt"
    path.write_text("".join(f"{5000 + 0.01 * i:.3f} 1.0\n" for i in range(400)))
    from tumult.spectrum_io import load_file
    return ctx.canvas.add_object(load_file(path), name="obs", source=("file", str(path)))


# --- the headless refresh -----------------------------------------------------------------
def test_refresh_ews_forwards_arguments_assigns_widths_and_names_failures(monkeypatch):
    calls = _stub_measure(monkeypatch, fail_at=(5002.0,))
    spec = Spectrum(np.arange(5000., 5010., 0.01), np.ones(1000))
    lines = [_line(5000.5, 5001.5, include=False, flags=("blended",)),
             _line(5002.0, 5003.0),
             _line(float("nan"), float("nan"))]
    before = [(ln.ew, ln.ew_err) for ln in lines]

    outcome = ab.refresh_ews(spec, lines, profile="voigt", snr=80.0, cont_pct=1.5)

    assert outcome.changed == (0,) and outcome.skipped == (2,)
    assert list(outcome.failed) == [1] and "no profile" in outcome.failed[1]
    assert [c["window"] for c in calls] == [(5000.5, 5001.5), (5002.0, 5003.0)]
    assert all(c["spec"] is spec and c["profile"] == "voigt" and c["snr"] == 80.0
               and c["cont_pct"] == 1.5 and c["x0"] == lines[i].lam0
               for i, c in enumerate(calls))
    assert (lines[0].ew, lines[0].ew_err) == (5000.5, 80.0)
    assert [(ln.ew, ln.ew_err) for ln in lines[1:]] == before[1:]   # failed and skipped keep
    assert lines[0].include is False and lines[0].flags == ("blended",)
    assert lines[0].fit_window == (5000.5, 5001.5)


def test_abundance_refresh_uses_the_stored_uncertainties_and_stops_after_synthesis(
        tmp_path, monkeypatch):
    from tumult.abundance import AbundanceResult
    calls = _stub_measure(monkeypatch, fail_at=(5002.0,))
    document, ctx = _build_document()
    messages = []
    ctx.notify = lambda text, level="info", **_: messages.append(text)
    obs = _observation(ctx, tmp_path)
    abund = ctx.tools._tools["abundance"]
    res = AbundanceResult(element="Fe", symbol="Fe", summary={}, uid="res-a", name="Fe A",
                          ref_uid=obs.uid, ref_name="obs", stage="lines", profile="gauss",
                          snr=120.0, cont_pct=0.7,
                          lines=[_line(5000.5, 5001.5, include=False, flags=("weak",)),
                                 _line(5002.0, 5003.0)])
    abund.results.append(res)
    abund.form1.set_values(snr=999.0, cont_pct=9.0)    # the live form is not the source

    ctx.emit("filter-state", exists=True, visible=True)
    assert not calls                                   # an inactive tool does not measure

    abund.select_result(res.uid)
    ctx.emit("filter-state", exists=True, visible=True)
    assert [c["window"] for c in calls] == [(5000.5, 5001.5), (5002.0, 5003.0)]
    assert all(c["snr"] == 120.0 and c["cont_pct"] == 0.7 and c["profile"] == "gauss"
               for c in calls)
    assert (res.lines[0].ew, res.lines[0].ew_err) == (5000.5, 120.0)
    assert (res.lines[1].ew, res.lines[1].ew_err) == (50.0, 2.0)      # the failed line
    assert res.lines[0].include is False and res.lines[0].flags == ("weak",)
    assert messages and "5002.50" in messages[-1] and "1 line" in messages[-1]

    calls.clear()
    for stage in ("synth", "uncert"):
        res.stage = stage
        ctx.emit("filter-state", exists=True, visible=True)
    assert not calls                                   # inverted abundances stand as measured


def test_balance_refresh_uses_the_stored_uncertainties_and_stops_once_solved(
        tmp_path, monkeypatch):
    from tumult import measurement as ms
    calls = _stub_measure(monkeypatch)
    document, ctx = _build_document()
    obs = _observation(ctx, tmp_path)
    ionbal = ctx.tools._tools["ionbal"]
    res = ms.BalanceResult(element="Fe", uid="res-b", name="Fe balance", ref_uid=obs.uid,
                           ref_name="obs", profile="voigt", snr=45.0, cont_pct=2.0,
                           lines=[_line(5000.5, 5001.5), _line(5002.0, 5003.0, include=False)])
    ionbal.results.append(res)
    ionbal.form1.set_values(snr=999.0, cont_pct=9.0)   # the live form is not the source
    ionbal.select_result(res.uid)

    ctx.emit("filter-state", exists=True, visible=True)
    assert len(calls) == 2
    assert all(c["snr"] == 45.0 and c["cont_pct"] == 2.0 and c["profile"] == "voigt"
               for c in calls)
    assert [(ln.ew, ln.ew_err) for ln in res.lines] == [(5000.5, 45.0), (5002.0, 45.0)]
    assert res.lines[1].include is False

    calls.clear()
    res.stage = ms.STAGE_SOLVED
    ctx.emit("filter-state", exists=True, visible=True)
    assert not calls


# --- the browser half ---------------------------------------------------------------------
def _playwright_node_path():
    """Where node resolves ``playwright`` from: ``TUMULT_NODE_PATH`` or the npx cache."""
    explicit = os.environ.get("TUMULT_NODE_PATH")
    candidates = [Path(explicit)] if explicit else []
    candidates += sorted(Path.home().glob(".npm/_npx/*/node_modules"))
    for base in candidates:
        if (base / "playwright" / "package.json").exists():
            return base
    return None


HARNESS = r"""
const { chromium } = require('playwright');
const input = JSON.parse(require('fs').readFileSync(0, 'utf8'));
(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  page.on('pageerror', (e) => { console.error('pageerror', e); process.exitCode = 2; });
  await page.setContent('<!doctype html><html><body>' + input.body
    + '<script>window.sent = []; window.Bokeh = { documents: [{ get_model_by_name: function (name) {'
    + '  return { setv: function (v) { window.sent.push({ channel: name, value: v.value }); } };'
    + '} }] };</script>'
    + '<script>' + input.scripts.join('\n') + '</script></body></html>');
  const observed = [];
  for (const step of input.steps) {
    if (step.op === 'click') await page.click(step.selector);
    else if (step.op === 'select') await page.selectOption(step.selector, step.value);
    else if (step.op === 'html') await page.$eval(step.selector, (el, h) => { el.innerHTML = h; }, step.html);
    else if (step.op === 'mark') await page.evaluate((m) => window.sent.push({ mark: m }), step.mark);
    else if (step.op === 'observe') observed.push(await page.$$eval(step.selector,
        (els, attr) => els.map((el) => el.getAttribute(attr)), step.attr));
  }
  const sent = await page.evaluate(() => window.sent);
  console.log(JSON.stringify({ sent, observed }));
  await browser.close();
})().catch((e) => { console.error(e); process.exit(1); });
"""


def _run_browser(tmp_path, body, scripts, steps):
    node = shutil.which("node")
    base = _playwright_node_path()
    if node is None or base is None:
        pytest.skip("node with a playwright package (npx cache or TUMULT_NODE_PATH) is "
                    "required for the browser half")
    harness = tmp_path / "harness.js"
    harness.write_text(HARNESS)
    proc = subprocess.run([node, str(harness)], input=json.dumps(
        dict(body=body, scripts=scripts, steps=steps)), capture_output=True, text=True,
        env={**os.environ, "NODE_PATH": str(base)}, timeout=120)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def _step_tabs(prefix):
    """The template's own step tabs for ``prefix``, verbatim."""
    match = re.search(rf'<div class="step-tabs" data-tool="{prefix}"[^>]*>.*?</div>',
                      TEMPLATE.read_text(), re.S)
    assert match, prefix
    return match.group(0)


def _channel(document, name):
    model = document.select_one({"name": name})
    assert model is not None, name
    return model


def _send(document, sent, index):
    """Apply the ``index``-th recorded payload to its channel, as the page would."""
    record = sent[index]
    _channel(document, record["channel"]).value = record["value"]
    return record


def test_browser_events_route_through_the_shared_binder_to_each_tool(tmp_path, monkeypatch):
    from tumult.abundance import AbundanceResult
    from tumult import measurement as ms
    import store
    from abundwizard import ANALYSIS_SCRIPT
    from synthsetup import ATMOS_SCRIPT, INDIVIDUAL_GRID

    spaced = "my grid: 3D"
    store.set_section("synthesis", {
        "atmos_grids": [{"id": "g1", "name": "Test grid", "folder": "", "pattern": "",
                         "format": "Marcs"},
                        {"id": spaced, "name": "Spaced: name", "folder": "", "pattern": "",
                         "format": "Marcs"}],
        "atmos_grid_id": "g1", "ref_grid_id": "g1"})
    store.set_section("abundcorr_run", {"atmos_grid_id": INDIVIDUAL_GRID})
    document, ctx = _build_document()
    obs = _observation(ctx, tmp_path)
    abund, corr, ionbal = (ctx.tools._tools[k] for k in ("abundance", "abundcorr", "ionbal"))
    res_a = AbundanceResult(element="Fe", symbol="Fe", summary={}, uid="res-a",
                            name="Fe: A", ref_uid=obs.uid, ref_name="obs", stage="lines",
                            lines=[_line(5000.5, 5001.5), _line(5002.0, 5003.0)])
    res_a2 = AbundanceResult(element="Ti", symbol="Ti", summary={}, uid="res-a2",
                             name="Ti A", ref_uid=obs.uid, ref_name="obs", stage="lines",
                             lines=[_line(5001.0, 5001.2)])
    res_c = AbundanceResult(element="Ba", symbol="Ba", summary={}, uid="res-c", name="Ba corr",
                            atom_path="atom.ba", lines=[_line(5000.5, 5001.5)])
    res_b = ms.BalanceResult(element="Fe", uid="res-b", name="Fe balance", ref_uid=obs.uid,
                             ref_name="obs", lines=[_line(5000.5, 5001.5), _line(5002., 5003.)])
    abund.results += [res_a, res_a2]
    corr.results.append(res_c)
    ionbal.results.append(res_b)
    for tool, res in ((abund, res_a), (corr, res_c), (ionbal, res_b)):
        tool.selected_uid = res.uid
    synth_setup, abund_setup, corr_setup = ctx.synth_modal.setup, abund.setup, corr.setup

    body = "".join([
        f'<div id="abund-results">{abund.results_html()}</div>',
        f'<div id="atom-objects">{"".join(corr.result_rows(corr.results))}</div>',
        f'<div id="ionbal-results">{ionbal.results_html()}</div>',
        f'<div id="abund-table">{abund.table_html(res_a)}</div>',
        f'<div id="corr-table">{corr.table_html(res_c)}</div>',
        f'<div id="ionbal-table">{ionbal.table_html(res_b)}</div>',
        _step_tabs("abund"), _step_tabs("corr"), _step_tabs("ionbal"),
        f'<div id="synth-grid">{synth_setup.grid_ui_html()}</div>',
        f'<div id="abund-grid">{abund_setup.grid_ui_html()}</div>',
        f'<div id="corr-grid">{corr_setup.grid_ui_html()}</div>',
    ])
    assert "line-include-chk" in body and 'data-atmos-key="format"' in body
    replaced_table = ionbal.table_html(res_b).replace('data-line="1"', 'data-line="0"')
    steps = [
        {"op": "click", "selector": '[data-tool="abund"] .step-tab[data-step="2"]'},
        {"op": "click", "selector": '[data-tool="abund"] .step-tab[data-step="2"]'},
        {"op": "click", "selector": '[data-tool="corr"] .step-tab[data-step="2"]'},
        {"op": "click", "selector": '#abund-table tr[data-line="1"] .line-include-chk'},
        {"op": "observe", "selector": "#abund-table tr[data-line]", "attr": "class"},
        {"op": "click", "selector": '#corr-table tr[data-line="0"] td:last-child'},
        {"op": "observe", "selector": "#corr-table tr[data-line]", "attr": "class"},
        {"op": "click", "selector": '#ionbal-table .line-sort[data-sortcol="ew"]'},
        {"op": "click", "selector": '#ionbal-table .analysis-mask-chk'},
        {"op": "click", "selector": '#abund-results [data-result="res-a2"] .object-name'},
        {"op": "click", "selector": '#atom-objects [data-result="res-c"] .object-name'},
        {"op": "click", "selector": '#abund-results [data-result-remove="res-a"]'},
        {"op": "click", "selector": '#ionbal-results [data-result-remove="res-b"]'},
        {"op": "select", "selector": '#synth-grid select[data-atmos-key="grid"]', "value": spaced},
        {"op": "select", "selector": '#abund-grid select[data-atmos-key="mode"]', "value": "pick"},
        {"op": "select", "selector": '#corr-grid select[data-atmos-key="format"]', "value": "Text"},
        {"op": "mark", "mark": "replaced"},
        {"op": "html", "selector": "#ionbal-table", "html": replaced_table},
        {"op": "click", "selector": '#ionbal-table tr[data-line="0"] td:last-child'},
    ]
    out = _run_browser(tmp_path, body, [ATMOS_SCRIPT, ANALYSIS_SCRIPT], steps)
    sent = out["sent"]
    # analysis channels are colon-delimited, the atmosphere channels space-delimited
    kinds = [(r["channel"], r["value"].split(" " if "grid" in r["channel"] else ":")[0])
             for r in sent if "value" in r]
    assert kinds == [
        ("abund-actions", "step"), ("abund-actions", "step"), ("corr-actions", "step"),
        ("abund-actions", "toggle-line"), ("corr-actions", "jump-line"),
        ("ionbal-actions", "sort"), ("ionbal-actions", "mask-toggle"),
        ("abund-actions", "select"), ("corr-actions", "select"),
        ("abund-actions", "remove"), ("ionbal-actions", "remove"),
        ("synth-grid-actions", "grid"), ("abund-grid-actions", "mode"),
        ("corr-grid-actions", "format"), ("ionbal-actions", "jump-line")]
    values = [r["value"] for r in sent if "value" in r]
    assert values[0] != values[1]                        # a repeated click is a new payload
    assert len({v.split(":")[2] for v in values[:3]}) == 3
    abund_classes, corr_classes = out["observed"]
    assert not any("row-selected" in c for c in abund_classes)   # the checkbox jumped nowhere
    assert [("row-selected" in c) for c in corr_classes] == [True]
    assert values[-1].startswith("jump-line:0:")               # from the replaced table

    # --- the same payloads through the document's channels ----------------------------
    jumps, sorts = [], []
    monkeypatch.setattr(corr, "_jump_line", lambda token: jumps.append(("corr", token)))
    monkeypatch.setattr(ionbal, "_jump_line", lambda token: jumps.append(("ionbal", token)))
    values = iter(v for v in sent if "value" in v)

    def apply(kind):
        record = next(values)
        assert record["value"].startswith(kind), (kind, record)
        _channel(document, record["channel"]).value = record["value"]

    apply("step"); assert abund.step == 2
    abund.step = 1
    apply("step"); assert abund.step == 2                # the repeat lands too
    apply("step"); assert corr.step == 2
    assert res_a.lines[1].include
    apply("toggle-line"); assert res_a.lines[1].include is False and res_a.lines[0].include
    apply("jump-line"); assert jumps == [("corr", "0")]
    apply("sort"); assert ionbal.sort.col == "ew"
    assert ionbal.show_toggled_only is False
    apply("mask-toggle"); assert ionbal.show_toggled_only is True
    apply("select"); assert abund.selected_uid == "res-a2" and ctx.focus_kind == abund.kind
    apply("select"); assert corr.selected_uid == "res-c" and ctx.focus_kind == corr.kind
    apply("remove"); assert [r.uid for r in abund.results] == ["res-a2"]
    assert abund.selected_uid == "res-a2"
    apply("remove"); assert ionbal.results == [] and ionbal.selected_uid is None
    from htmlbridge import parse_action
    grid_payload = next(v for v in sent if v.get("channel") == "synth-grid-actions")["value"]
    assert parse_action(grid_payload) == ("grid", spaced)
    # atmosphere edits are a calculation draft: the setup reflects them at once, the
    # preference store only once the calculation is submitted (setup.commit)
    apply("grid")
    assert synth_setup.selected_grid_id() == spaced
    assert store.get_section("synthesis").get("atmos_grid_id") != spaced
    synth_setup.commit()
    assert store.get_section("synthesis")["atmos_grid_id"] == spaced
    apply("mode"); assert abund_setup.show_file()
    assert store.get_section("abund_run").get("atmos_mode") != "pick"
    abund_setup.commit(); assert store.get_section("abund_run")["atmos_mode"] == "pick"
    apply("format")
    assert "atmos_individual_format" not in store.get_section("abundcorr_run")
    corr_setup.commit()
    assert store.get_section("abundcorr_run")["atmos_individual_format"] == "Text"
    apply("jump-line"); assert jumps[-1] == ("ionbal", "0")
    assert next(values, None) is None


def test_analysis_and_atmosphere_scripts_are_emitted_once():
    """The page carries each shared script exactly once and no per-tool copies."""
    text = TEMPLATE.read_text()
    assert text.count("{{ analysis_script|safe }}") == 1
    assert text.count("{{ atmos_script|safe }}") == 1
    for stale in ("abundGridAction", "corrGridAction", "synthGridAction", "-grid-select",
                  "selectAnalysisRow", "data-corr-remove", "data-ionbal-remove"):
        assert stale not in text, stale
    document, ctx = _build_document()
    from abundwizard import ANALYSIS_SCRIPT
    from synthsetup import ATMOS_SCRIPT
    assert document.template_variables["analysis_script"] is ANALYSIS_SCRIPT
    assert document.template_variables["atmos_script"] is ATMOS_SCRIPT
