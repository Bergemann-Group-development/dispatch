"""The saved-workspace schema: validation, headless preparation and the partial-load
policy (TumultStudio2/session.py). Serialization only; no scientific content is checked."""
import copy
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "TumultStudio2"))

from tumult import Spectrum  # noqa: E402
from tumult.abundance import AbundanceResult, LineResult  # noqa: E402
from tumult.linemask import LineMask, LineMaskSet  # noqa: E402
from tumult.measurement import BalanceResult  # noqa: E402
import session  # noqa: E402
from session import SessionError, Workspace  # noqa: E402


def _line(lo=5000.5, hi=5001.5):
    mask = LineMask(lam0=0.5 * (lo + hi))
    mask.set_fit(lo, hi, 0.5)
    return LineResult(lam0=0.5 * (lo + hi), ion="Fe I", e_low=1.0, loggf=-1.0, ew=50.0,
                      mask=mask)


def _abundance(uid, ref_uid, atom_path=""):
    return AbundanceResult(element="Fe", symbol="Fe", summary={}, uid=uid, name=f"A {uid}",
                           ref_uid=ref_uid, ref_name="obs", lines=[_line()],
                           atom_path=atom_path)


def _balance(uid, ref_uid, maskset_uid):
    return BalanceResult(element="Fe", uid=uid, name=f"B {uid}", ref_uid=ref_uid,
                         ref_name="obs", maskset_uid=maskset_uid, maskset_name="masks",
                         lines=[_line()])


def _snapshot(tmp_path):
    """A complete, valid version-2 slot: two spectra, two mask sets, one atom, one
    result of each kind, a focused correction result."""
    obs = tmp_path / "obs.txt"
    obs.write_text("".join(f"{5000 + 0.01 * i:.3f} 1.0\n" for i in range(300)))
    synth = Spectrum(np.linspace(5000., 5003., 30), np.ones(30))
    edited = LineMaskSet(masks=[LineMask(lam0=5001., lmin=5000.5, lmax=5001.5,
                                         smin=5000., smax=5002.)],
                         pad=0.5, label="fe-lines", uid=2, element="Fe",
                         source_path=str(tmp_path / "fe-lines.txt"))
    edited.mark_clean()
    edited.masks[0].lmax = 5001.6          # edited away from its file
    atom = str(tmp_path / "atom.fe")
    return {
        "schema_version": session.SCHEMA_VERSION,
        "spectra": [
            {"uid": "s-file", "kind": "file", "name": "obs", "path": str(obs),
             "color": "#ffffff", "dash": "solid", "width": 1.5, "alpha": 0.9,
             "visible": False, "conv": {"vshift": 1.0},
             "norm_knots": {"x": [5000.0, 5002.0], "y": [1.0, 1.0]}},
            {"uid": "s-synth", "kind": "synth", "name": "model", "data": synth.to_dict(),
             "meta": {"teff": 5000}, "already_normalised": True}],
        "spectrum_selected": "s-file",
        "masks": {"sets": [LineMaskSet(uid=1).to_dict(), edited.to_dict()],
                  "current": 2, "active": 0},
        "atoms": [{"uid": "a-1", "path": atom}],
        "atom_selected": "a-1",
        "results": {
            "abundance": {"results": [_abundance("r-a", "s-file").to_dict()],
                          "selected": "r-a"},
            "abundcorr": {"results": [_abundance("r-c", None, atom).to_dict()],
                          "selected": "r-c"},
            "ionbal": {"results": [_balance("r-b", "s-file", "2").to_dict()],
                       "selected": None}},
        "view": {"xunit": "Å", "yunit": "normalized", "main_x": [5000.2, 5001.8],
                 "main_y": [0.4, 1.1], "nav_x": [5000.0, 5003.0], "nav_y": None},
        "filter": {"active": True, "method": "lowpass", "res": 20000.0},
        "focus": {"kind": "corr", "uid": "r-c"},
        "tools": {"active": "abundcorr", "atom_tab": "trans",
                  "steps": {"abundcorr": 2}, "jumped": {"abundcorr": 0}},
    }


def _tools():
    """The wizards' persistence hooks without a document (class-level factories)."""
    from abundancetool import AbundanceTool
    from abundcorrtool import AbundCorrTool
    from ionbalance import IonBalanceTool
    return {name: cls.__new__(cls) for name, cls in (
        ("abundance", AbundanceTool), ("abundcorr", AbundCorrTool), ("ionbal", IonBalanceTool))}


def _prepare(ws, load_atom=lambda path: object()):
    from tumult.spectrum_io import load_file

    def load_spectrum(snap):
        return Spectrum.from_dict(snap.data) if snap.kind == "synth" else load_file(snap.path)
    return session.prepare(ws, load_spectrum=load_spectrum, load_atom=load_atom,
                           tools=_tools())


# --- schema ---------------------------------------------------------------------------
def test_workspace_round_trips_through_its_dict(tmp_path):
    snap = _snapshot(tmp_path)
    ws = Workspace.from_dict(snap)
    again = Workspace.from_dict(ws.to_dict())
    assert again == ws
    assert ws.to_dict()["spectra"][0]["norm_knots"] == {"x": [5000.0, 5002.0], "y": [1.0, 1.0]}
    assert ws.focus.kind == "corr" and ws.focus.uid == "r-c"
    assert ws.uids_of("linemask") == [1, 2] and ws.uids_of("corr") == ["r-c"]


@pytest.mark.parametrize("version", [1, None, "2", 3])
def test_unsupported_versions_are_rejected_before_anything_else(tmp_path, version):
    snap = _snapshot(tmp_path)
    snap["schema_version"] = version
    snap["spectra"] = "not even a list"        # never reached: the version comes first
    with pytest.raises(SessionError, match=f"version {version!r} is not supported"):
        Workspace.from_dict(snap)
    assert "unsupported format" in session.summary(snap)


def _broken(tmp_path):
    base = _snapshot(tmp_path)

    def variant(mutate):
        snap = copy.deepcopy(base)
        mutate(snap)
        return snap

    def bad_knots(s):
        s["spectra"][0]["norm_knots"] = {"x": [1.0, 2.0], "y": [1.0]}

    def nan_knots(s):
        s["spectra"][0]["norm_knots"] = {"x": [1.0, float("nan")], "y": [1.0, 1.0]}

    def dup(s):
        s["spectra"][1]["uid"] = "s-file"

    def no_uid(s):
        del s["spectra"][0]["uid"]

    def bad_kind(s):
        s["spectra"][0]["kind"] = "web"

    def synth_no_data(s):
        del s["spectra"][1]["data"]

    def bool_width(s):
        s["spectra"][0]["width"] = True

    def focus_unknown(s):
        s["focus"] = {"kind": "spectrum", "uid": "gone"}

    def focus_kind(s):
        s["focus"] = {"kind": "planet", "uid": None}

    def selected_unknown(s):
        s["results"]["abundance"]["selected"] = "nope"

    def mask_current(s):
        s["masks"]["current"] = 9

    def mask_uid(s):
        s["masks"]["sets"][0]["uid"] = 0

    def result_no_uid(s):
        del s["results"]["ionbal"]["results"][0]["uid"]

    def bad_step(s):
        s["tools"]["steps"] = {"abundance": True}

    def bad_range(s):
        s["view"]["main_x"] = [5001.0, 5000.0]

    def atom_selected(s):
        s["atom_selected"] = "a-9"

    def unknown_section(s):
        s["results"]["fits"] = {}

    def wrong_shape(s):
        s["masks"] = []

    return {fn.__name__: (variant(fn), where) for fn, where in (
        (bad_knots, "norm_knots"), (nan_knots, "norm_knots"), (dup, "duplicate spectrum uid"),
        (no_uid, "missing 'uid'"), (bad_kind, "unknown kind"), (synth_no_data, "embedded data"),
        (bool_width, "width"), (focus_unknown, "focus"), (focus_kind, "focus kind"),
        (selected_unknown, "selected result"), (mask_current, "current set"),
        (mask_uid, "must be positive"), (result_no_uid, "ionbal"), (bad_step, "steps"),
        (bad_range, "main_x"), (atom_selected, "atom_selected"),
        (unknown_section, "unknown section"), (wrong_shape, "masks"))}


def test_malformed_slots_name_the_offending_field(tmp_path):
    for name, (snap, where) in _broken(tmp_path).items():
        with pytest.raises(SessionError, match=where):
            Workspace.from_dict(snap)


def test_summary_describes_a_slot_without_loading_it(tmp_path):
    snap = _snapshot(tmp_path)
    assert session.summary(snap) == "2 spectra · 1 line mask · 1 atom · 3 results"
    assert session.summary({"schema_version": 2, "spectra": []}) == "0 spectra"
    assert session.summary(None) == "unreadable"


def test_scientific_records_accept_only_their_own_version():
    from tumult.measured import MeasuredLines
    record = _abundance("r", "s").to_dict()
    assert AbundanceResult.from_dict(record).uid == "r"
    for bad in (dict(record, schema_version=0), {k: v for k, v in record.items()
                                                  if k != "schema_version"}, "text"):
        with pytest.raises(ValueError, match="schema version|must be an object"):
            MeasuredLines.from_dict(bad)


# --- preparation and the partial-load policy -----------------------------------------
def test_prepare_builds_every_collection_off_the_workspace(tmp_path):
    ws = Workspace.from_dict(_snapshot(tmp_path))
    atoms = []
    prepared = _prepare(ws, load_atom=lambda path: atoms.append(path) or {"path": path})
    assert not prepared.skipped and not prepared.unresolved
    assert [p.snap.uid for p in prepared.spectra] == ["s-file", "s-synth"]
    assert prepared.spectra[1].meta == {"teff": 5000}
    assert [s.uid for s in prepared.mask_sets] == [1, 2]
    edited = prepared.mask_sets[1]
    assert edited.dirty and edited.masks[0].lmax == 5001.6      # edit and baseline survive
    assert edited.revert() and edited.masks[0].lmax == 5001.5
    assert atoms == [str(tmp_path / "atom.fe")] and prepared.atom_selected == "a-1"
    assert {k: (len(r), sel) for k, (r, sel) in prepared.results.items()} == {
        "abundance": (1, "r-a"), "abundcorr": (1, "r-c"), "ionbal": (1, None)}
    assert prepared.results["ionbal"][0][0].maskset_uid == "2"
    assert prepared.focus.kind == "corr" and prepared.focus.uid == "r-c"
    assert prepared.counts() == {"spectra": 2, "masks": 1, "atoms": 1, "results": 3}


def test_missing_files_are_skipped_with_reasons_and_their_links_cleared(tmp_path):
    snap = _snapshot(tmp_path)
    snap["spectra"][0]["path"] = str(tmp_path / "gone.txt")
    snap["focus"] = {"kind": "spectrum", "uid": "s-file"}

    def no_atom(path):
        raise FileNotFoundError(path)
    prepared = _prepare(Workspace.from_dict(snap), load_atom=no_atom)
    assert [(s.kind, s.uid, s.path) for s in prepared.skipped] == [
        ("spectrum", "s-file", str(tmp_path / "gone.txt")),
        ("atom", "a-1", str(tmp_path / "atom.fe"))]
    assert all("FileNotFoundError" in s.reason for s in prepared.skipped)
    # the independent objects are restored, the references to the missing ones cleared
    assert [p.snap.uid for p in prepared.spectra] == ["s-synth"]
    assert len(prepared.mask_sets) == 2 and not prepared.atoms
    abund, balance = prepared.results["abundance"][0][0], prepared.results["ionbal"][0][0]
    assert abund.ref_uid is None and balance.ref_uid is None and balance.maskset_uid == "2"
    assert len(prepared.unresolved) == 3 and any("atom.fe" in u for u in prepared.unresolved)
    # nothing is selected in place of the missing objects
    assert prepared.spectrum_selected is None and prepared.atom_selected is None
    assert prepared.focus.kind == "spectrum" and prepared.focus.uid is None


def test_a_missing_mask_set_clears_the_balance_link(tmp_path):
    snap = _snapshot(tmp_path)
    snap["masks"] = {"sets": [], "current": None, "active": None}
    prepared = _prepare(Workspace.from_dict(snap))
    balance = prepared.results["ionbal"][0][0]
    assert balance.maskset_uid == "" and any("line mask" in u for u in prepared.unresolved)


def test_broken_embedded_objects_abort_preparation(tmp_path):
    good = _snapshot(tmp_path)
    for mutate, where in (
            (lambda s: s["results"]["abundance"]["results"][0].__setitem__("lines", [3]),
             "results.abundance"),
            (lambda s: s["results"]["ionbal"]["results"][0].__setitem__("schema_version", 0),
             "results.ionbal"),
            (lambda s: s["masks"]["sets"][1].__setitem__("pad", "wide"), r"masks.sets\[1\]"),
            (lambda s: s["spectra"][1].__setitem__("data", {"x": "no"}), "embedded data")):
        snap = copy.deepcopy(good)
        mutate(snap)
        with pytest.raises(SessionError, match=where):
            _prepare(Workspace.from_dict(snap))


def test_empty_and_single_collection_slots_prepare(tmp_path):
    empty = {"schema_version": session.SCHEMA_VERSION}
    prepared = _prepare(Workspace.from_dict(empty))
    assert prepared.counts() == {"spectra": 0, "masks": 0, "atoms": 0, "results": 0}
    assert prepared.focus.kind == "spectrum" and prepared.focus.uid is None
    full = _snapshot(tmp_path)
    only_atom = {"schema_version": 2, "atoms": full["atoms"], "atom_selected": "a-1",
                 "focus": {"kind": "atom", "uid": "a-1"}}
    prepared = _prepare(Workspace.from_dict(only_atom), load_atom=lambda p: p)
    assert prepared.counts()["atoms"] == 1 and prepared.focus.uid == "a-1"
    only_results = {"schema_version": 2, "results": full["results"],
                    "focus": {"kind": "ionbal", "uid": "r-b"}}
    prepared = _prepare(Workspace.from_dict(only_results))
    assert prepared.counts()["results"] == 3 and prepared.focus.uid == "r-b"
    assert len(prepared.unresolved) == 4      # two observations, one mask set, one atom
