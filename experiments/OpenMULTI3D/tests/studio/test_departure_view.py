"""Solving and looking at departure coefficients from a model atom.

What matters to someone using it: the atom's panel starts one solve of one atmosphere,
which is listed under that atom the moment it is submitted and stays there while it runs;
an interpolated grid is refused, because that would be a departure grid and not this
object; and once the populations exist, the row selected in the table below the plot is
what is drawn — a level with its departure relative to its own stage's ground state, a
transition with both its levels and their ratio.

The solve itself is not run here: the submission is captured instead of reaching the
workflow daemon, and the populations are a stub of
:meth:`tumult.departures.Departures.profile`, whose own reading is covered in
``tests/io/test_departure_profile.py``. The model atom is real, so the level and
transition indices, ionisation stages and ground states are the file's own.
"""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("bokeh")

ROOT = Path(__file__).resolve().parents[2]
ATOM = ROOT / "input" / "atoms" / "atom.ba06"
pytestmark = pytest.mark.skipif(not ATOM.is_file(),
                                reason="input/atoms/atom.ba06 is absent")


class _Submission:
    """Stands in for ``ctx.submit_operation``: keeps the callbacks a solve hands over."""

    def __init__(self):
        self.calls = []

    def __call__(self, name, operation, request, on_done, *, kwargs=None, meta=None,
                 on_progress=None, on_error=None, on_native=None):
        from tumult.synthesis.cancel import CancelToken

        class Job:
            pid = len(self.calls) + 1
            token = CancelToken()

        self.calls.append(dict(name=name, operation=operation, request=request,
                               kwargs=kwargs, done=on_done, native=on_native,
                               error=on_error))
        return Job()


def _request():
    from tumult.synthesis import Atmosphere, Star, SynthRequest
    from tumult.synthesis.inputs import SynthConfig
    return SynthRequest(config=SynthConfig(root=ROOT),
                        atmosphere=Atmosphere("input/atmos/marcs.mod", "Marcs"),
                        star=Star(5777.0, 4.44, 0.0, 1.0))


def _solved(tmp_path, atom, profile):
    """A solved run of ``atom``, with its populations already read."""
    from tumult.departures import Departures
    departures = Departures(folder=tmp_path, name=Path(tmp_path).name,
                            atmosphere="marcs.mod", species=("Ba",), atom=str(atom),
                            star={"teff": 5777.0, "logg": 4.44, "feh": 0.0, "vmic": 1.0})
    object.__setattr__(departures, "profile", lambda species=None: profile)
    return departures


@pytest.fixture
def studio(bokeh_document, studio_config):
    from TumultStudio import app
    return app.build()


@pytest.fixture
def atom(studio):
    from tumult.atom import Atom
    return studio.tools.get("atom").add(str(ATOM), Atom.read(ATOM))


@pytest.fixture
def profile(atom):
    """A departure profile of the real atom's shape: every level in LTE except level 3,
    which is depleted towards the surface, on three depth points."""
    from tumult.departures import DepartureProfile
    b = np.ones((atom.atom.nk, 3))
    b[3] = [0.25, 0.5, 1.0]
    return DepartureProfile("Ba", np.array([-2.0, 0.0, 2.0]), b, columns=1)


def test_the_atom_panel_starts_one_solve_and_lists_it_while_it_runs(
        studio, atom, monkeypatch, tmp_path, profile):
    ctx = studio
    tool, atomtool = ctx.tools.get("departures"), ctx.tools.get("atom")
    submission = _Submission()
    monkeypatch.setattr(ctx, "submit_operation", submission)
    monkeypatch.setattr(ctx, "schedule", lambda fn: fn())     # the deferred submission

    assert "atom-departures" in atomtool.panel_html()         # the entry point exists
    tool.open_for(atom)
    assert "Ba" in tool.context_html()
    monkeypatch.setattr(tool.setup, "context", _request)
    tool.run()

    # the object is listed under its atom, selected and in view before anything is computed
    [obj] = tool.objects
    assert tool.selected() is obj and ctx.workspace.focus_kind == tool.kind
    assert obj.name.startswith("Ba departures") and tool.is_computing(obj)
    rows = atomtool.objects_html()
    assert f'data-result="{obj.uid}"' in rows and "Solving Ba departures" in rows
    # the solve carries the atom it was opened from, and nothing else
    [call] = submission.calls
    assert call["operation"] == "departures.solve"
    assert call["kwargs"] == {"atom": str(ATOM), "line_grid": True}
    assert call["request"].windows == () and call["request"].departures is None

    # how far the running dispatch run is reaches the sidebar row
    call["native"](0.5)
    assert "50%" in atomtool.objects_html()

    # the solved run lands: the object holds it and the strip is gone
    call["done"](_solved(tmp_path, ATOM, profile))
    assert obj.solved and not obj.progress
    assert "marcs.mod" in tool.panel_html() and "Use in synthesis" in tool.panel_html()


def test_the_panel_button_reaches_the_modal_through_the_page_script(studio, atom):
    """The atom panel's action, from the markup it renders to the modal it opens.

    The page routes a click on a panel button to the atom's action channel by its class
    (``static/js/analysis.js``), and the tool acts on the payload that arrives. A button
    the script does not route is inert in the browser while every server-side call still
    works, which is exactly what this checks.
    """
    atomtool, tool = studio.tools.get("atom"), studio.tools.get("departures")
    atomtool.select(atom.uid)
    assert "atom-departures" in atomtool.panel_html()
    script = (ROOT / "TumultStudio" / "static" / "js" / "analysis.js").read_text()
    assert ".atom-departures" in script and "atomAction('departures')" in script

    # what that click sends on the atom channel: "<kind> <arg> <nonce>", no argument here
    atomtool.actions.value = "departures  1700000000000"

    assert tool._pending_atom is atom
    assert tool._modal.value.startswith("open")
    assert "Ba" in tool._context_bridge.value


def test_an_interpolated_grid_is_refused_with_a_reason(studio, atom, monkeypatch):
    from tumult.synthesis import AtmosphereGrid
    ctx = studio
    tool = ctx.tools.get("departures")
    submission = _Submission()
    monkeypatch.setattr(ctx, "submit_operation", submission)
    monkeypatch.setattr(ctx, "schedule", lambda fn: fn())
    messages = []
    monkeypatch.setattr(ctx, "notify", lambda text, **kw: messages.append(text))

    tool.open_for(atom)
    grid = _request().update(atmosphere=AtmosphereGrid("input/atmos", "Marcs"))
    monkeypatch.setattr(tool.setup, "context", lambda: grid)
    tool.run()

    assert tool.objects == [] and submission.calls == []
    assert "one atmosphere" in messages[-1]


def test_the_modal_picks_one_model_and_never_interpolates_a_grid(studio, atom):
    """A departure solve is of one atmosphere, so the atmosphere cluster offers no
    interpolation: a registered grid is a list of models to pick from. The submit-time
    check stays for a request formed another way (`test_an_interpolated_grid_is_refused`),
    but nothing in the modal can ask for one."""
    from TumultStudio.synthesis import atmospicker
    setup = studio.tools.get("departures").setup
    setup.grid_ui_html()                      # a grid is listed to pick a model from
    assert 'data-atmos-key="grid"' in setup.grid_ui_html()
    # no Source toggle: the row states the only source there is
    assert 'data-atmos-key="mode"' not in setup.grid_ui_html()
    assert "Pick model" in setup.grid_ui_html()
    assert setup.show_file()                  # so the picker is open from the start

    # and it is the grid's parsed-model picker, not the path browser: a visible picker is
    # in node mode from the moment the scan starts, and says so while the grid is read
    picker = setup.atmos_picker
    assert picker.node_mode and picker.root is None
    assert "Filter models" in picker.query.placeholder
    assert picker._problem == atmospicker.SCANNING or picker._nodes

    # an interpolate choice arriving from the page anyway changes nothing
    setup.grid_actions.value = "mode interp 1700000000000"
    assert setup.show_file()

    # the synthesis modal's own cluster still interpolates
    assert 'data-atmos-key="mode"' in studio.synth_modal.setup.grid_ui_html()


def test_the_selected_level_and_transition_are_what_is_drawn(
        studio, atom, tmp_path, profile):
    from TumultStudio.atoms.departuretool import DepartureObject
    tool = studio.tools.get("departures")
    obj = DepartureObject(atom_path=str(ATOM), name="Ba departures", species="Ba",
                          departures=_solved(tmp_path, ATOM, profile))
    tool.add(obj)
    tool._profiles[obj.uid] = profile

    # a level: its own departure, and its departure relative to its stage's ground state
    tool.set_tab("levels")
    tool.pick_row("3")
    labels = list(tool.src.data["label"])
    assert len(labels) == 2 and labels[0].startswith("b — ")
    assert labels[1].startswith("b / b(ground Ba ")
    ground = tool._ground_of(tool.atom_of(obj), 3)
    assert list(tool.src.data["key"]) == ["b₄", f"b₄ / b{'₁₂₃₄₅₆₇₈₉'[ground]}"]
    assert list(tool.src.data["ys"][0]) == [0.25, 0.5, 1.0]
    assert list(tool.src.data["xs"][0]) == [-2.0, 0.0, 2.0]
    # the ground state of that stage is in LTE here, so the ratio is the level's own b
    assert list(tool.src.data["ys"][1]) == [0.25, 0.5, 1.0]
    # a ground state has no ratio of its own to show
    tool.pick_row("0")
    assert len(tool.src.data["label"]) == 1

    # a transition: both its levels and the ratio that is its line source function
    transition = next(i for i, tr in enumerate(atom.lines) if 3 in (tr.i, tr.j))
    tool.set_tab("trans")
    tool.pick_row(str(transition))
    labels = list(tool.src.data["label"])
    assert len(labels) == 3
    assert labels[0].startswith("b lower") and labels[1].startswith("b upper")
    assert labels[2].startswith("b_u / b_l")
    lower, upper = atom.lines[transition].i, atom.lines[transition].j
    assert np.allclose(tool.src.data["ys"][2],
                       profile.b[upper] / profile.b[lower])

    # the legend names each curve by level number as the tables count them; the full
    # text is what the hover shows
    from TumultStudio.atoms.departuretool import b_of
    b_l, b_u = b_of(lower), b_of(upper)
    assert b_of(11) == "b₁₂"
    keys = [b_l, b_u, f"{b_u} / {b_l}"]
    assert list(tool.src.data["key"]) == keys
    assert [item.label.value for item in tool.legend.items] == keys

    # a pin of the selected row adds nothing: the curves are drawn once
    tool.pin()
    assert len(tool.legend.items) == 3
    # once another row is selected the pin keeps its curves, each with its key in the
    # legend after the selection's own, in a run colour of its own
    tool.pick_row("0")                                  # another transition of the same atom
    assert tool.legend.visible and len(tool.legend.items) == 6
    assert [item.label.value for item in tool.legend.items[3:]] == keys
    assert not any("pinned" in item.label.value for item in tool.legend.items)
    # the selection's curves are drawn first, each on a renderer of its own with the dash
    # that tells the quantities apart; the pinned ones follow, in another colour, and say
    # so on hover
    drawn = [r for _s, r in tool._lines[:6]]
    assert drawn[0].glyph.line_color != drawn[3].glyph.line_color
    assert all(r.visible for r in drawn) and not any(r.visible for _s, r in tool._lines[6:])
    solid, dashed, dotted = [], [6], [2, 4]                  # Bokeh's named patterns
    assert [r.glyph.line_dash for r in drawn[:3]] == [solid, dashed, dotted]
    assert [r.glyph.line_dash for r in drawn[3:]] == [solid, dashed, dotted]
    assert all(source.data["label"][0].startswith("pinned · ")
               for source, _r in tool._lines[3:6])

    # clearing the pins takes the entry and its curves away
    tool.actions.value = "unpin::1700000000000"
    assert len(tool.legend.items) == 3
    assert all(not isinstance(k, tuple) for k in tool._colors)   # the pin's colour is freed
    assert not any(r.visible for _s, r in tool._lines[3:])


def test_the_table_below_the_plot_is_the_atom_s_own(studio, atom, tmp_path, profile):
    from TumultStudio.atoms.departuretool import DepartureObject
    tool = studio.tools.get("departures")
    obj = DepartureObject(atom_path=str(ATOM), name="Ba departures", species="Ba",
                          departures=_solved(tmp_path, ATOM, profile))
    tool.add(obj)
    tool.pick_row("3")
    tool.render_tables()

    levels = tool._tables["levels"].value
    assert 'data-line="3"' in levels and "row-selected" in levels
    assert 'data-line="0"' in tool._tables["trans"].value
    # the search narrows the table the selection is made in
    tool.query["levels"] = "no-such-level"
    tool.render_table("levels")
    assert "No matches" in tool._tables["levels"].value


def test_opening_a_solve_draws_the_first_row_of_the_table(studio, atom, tmp_path, profile):
    """A view that opened on no row would show an empty plot, so the first row of the
    table as displayed — its sort order and its search — is what is drawn."""
    from TumultStudio.atoms.atomtool import level_order
    from TumultStudio.atoms.departuretool import DepartureObject
    tool = studio.tools.get("departures")
    obj = DepartureObject(atom_path=str(ATOM), name="Ba departures", species="Ba",
                          departures=_solved(tmp_path, ATOM, profile))
    tool._profiles[obj.uid] = profile
    tool.add(obj)                                  # selects it and opens the view

    first = level_order(atom, tool.sort["levels"], "")[0]
    assert tool.row["levels"] == first
    assert tool.src.data["label"] and "row-selected" in tool._tables["levels"].value

    # sorting moves the table; a selection still shown is kept
    tool.sort["levels"].toggle("ecm")
    tool.render_table("levels")
    assert tool.row["levels"] == first

    # a search that hides it moves the selection to the first row that is shown
    shown = [i for i in level_order(atom, tool.sort["levels"], "II") if i != first]
    assert shown, "the test atom has no second ionisation stage to filter to"
    tool.actions.value = "search:levels II:1700000000000"
    assert tool.row["levels"] == level_order(atom, tool.sort["levels"], "II")[0]
    assert tool.row["levels"] != first
    assert tool.src.data["label"]                  # and the curves followed it


def test_a_solve_survives_the_session_as_its_run_folder(
        studio, atom, tmp_path, profile, monkeypatch):
    """The slot stores the run, not the populations: a restored solve reads them again."""
    from TumultStudio.atoms.departuretool import DepartureObject
    from TumultStudio.workspace import session
    from TumultStudio.workspace.session import Workspace

    tool = studio.tools.get("departures")
    departures = _solved(tmp_path, ATOM, profile)
    obj = DepartureObject(atom_path=str(ATOM), name="Ba departures", species="Ba",
                          departures=departures)
    tool.add(obj)

    saved = tool.snapshot()
    assert saved["departures"] == [{"uid": obj.uid, "name": "Ba departures",
                                    "atom_path": str(ATOM), "folder": str(tmp_path),
                                    "species": "Ba"}]
    assert saved["selected"] == obj.uid

    # the slot validates, and prepares by reopening the run
    ws = Workspace.from_dict({"schema_version": session.SCHEMA_VERSION,
                              "departures": saved["departures"],
                              "departure_selected": saved["selected"],
                              "focus": {"kind": "departures", "uid": obj.uid}})
    prepared = session.prepare(ws, load_spectrum=lambda snap: None,
                               load_atom=lambda path: None,
                               load_departures=lambda folder: departures,
                               deserialise=lambda name, records: [],
                               relink=lambda *a, **k: [])
    assert [p.snap.uid for p in prepared.departures] == [obj.uid]
    assert prepared.focus.kind == "departures"

    tool.restore(prepared.departures, prepared.departure_selected)
    [restored] = tool.objects
    assert restored.uid == obj.uid and restored.departures is departures
    assert restored.name == "Ba departures" and restored.atom_path == str(ATOM)

    # a run that is gone is skipped with its reason rather than failing the load
    def missing(folder):
        raise FileNotFoundError(f"No departures record or params.nml in {folder}")

    prepared = session.prepare(ws, load_spectrum=lambda snap: None,
                               load_atom=lambda path: None, load_departures=missing,
                               deserialise=lambda name, records: [],
                               relink=lambda *a, **k: [])
    assert prepared.departures == [] and prepared.focus.kind == "spectrum"
    [skipped] = prepared.skipped
    assert skipped.kind == "departures" and "params.nml" in skipped.reason


def test_removing_the_atom_takes_its_solves_and_undo_puts_them_back(
        studio, atom, tmp_path, profile):
    from TumultStudio.atoms.departuretool import DepartureObject
    ctx = studio
    tool, atomtool = ctx.tools.get("departures"), ctx.tools.get("atom")
    obj = DepartureObject(atom_path=str(ATOM), name="Ba departures", species="Ba",
                          departures=_solved(tmp_path, ATOM, profile))
    tool.add(obj)

    plan = ctx.removal.plan("atom", atom.uid)
    assert plan.dependents == ["Ba departures"] and plan.needs_confirmation

    ctx.removal.request("atom", atom.uid)
    ctx.removal.confirm(next(iter(ctx.removal.pending)))
    assert tool.objects == [] and atomtool.atoms == []

    ctx.removal.undo(ctx.removal._undo[0])
    assert [o.uid for o in tool.objects] == [obj.uid]
    assert [a.uid for a in atomtool.atoms] == [atom.uid]
