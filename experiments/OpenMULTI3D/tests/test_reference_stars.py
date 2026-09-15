"""The reference-star catalogue and the grid-model picker of the synthesis setups.

The catalogue is one list (tumult.refstars) three consumers read; a setup's dropdown and
its four parameter fields must agree in both directions; a registered grid is picked from
by parsed model, not by browsing its folder.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'TumultStudio2'))

import store
from tumult import linestrength, refstars

MARCS_NAMES = ["p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00.mod",
               "p5750_g+4.5_m0.0_t02_st_z+0.00_a+0.00.mod",
               "p4500_g+2.0_m1.0_t02_st_z-1.00_a+0.40.mod",
               "s3800_g+0.0_m15._t05_st_z+0.00_a+0.00.mod"]


@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(store, "CONFIG_PATH", tmp_path / "prefs.json")
    saved = refstars.stars()
    refstars.set_stars(refstars.DEFAULT_STARS)
    yield
    refstars.set_stars(saved)
    linestrength.set_marcs_grid(None)


def _marcs_folder(tmp_path):
    folder = tmp_path / "marcs"
    folder.mkdir()
    for name in MARCS_NAMES:
        (folder / name).write_text("")
    return folder


def _setup(tmp_path, folder):
    from gridmanager import GridManager
    from synthsetup import SynthSetup
    store.set_section("synthesis", {
        "atmos_grids": [{"id": "g1", "name": "Tiny MARCS", "folder": str(folder),
                         "pattern": "", "format": "Marcs"}],
        "atmos_grid_id": "g1"})
    store.set_section("test_run", {"atmos_mode": "pick", "atmos_grid_id": "g1"})
    grids = GridManager(default_pattern="", default_storage=tmp_path / "none")
    return SynthSetup("tst", "test_run", grids)


def test_dropdown_and_fields_stay_in_step(tmp_path):
    setup = _setup(tmp_path, _marcs_folder(tmp_path))
    changes = []
    setup.on_run_change = lambda k, v: changes.append(k)
    assert setup.run_form.value("refstar") == refstars.CUSTOM_KEY   # 5777/4.44 is nobody's

    setup.run_form._on_input("value", "", "refstar=rgb 1")
    assert setup.stellar_values() == refstars.get("rgb").params
    assert store.get_section("test_run") == {"atmos_mode": "pick", "atmos_grid_id": "g1"}
    assert set(changes) == set(refstars.PARAM_KEYS)             # the host saw every field

    setup.run_form._on_input("value", "", "teff=4600 2")
    assert setup.run_form.value("refstar") == refstars.CUSTOM_KEY
    assert setup.run_form.value("teff") == 4600

    # values written back (a fit result, a picked file) land on the matching entry
    setup.set_stellar({"teff": 4500, "logg": 2.0, "feh": 0.0, "vmic": 2.0})
    assert setup.run_form.value("refstar") == "rgb"
    setup.commit()
    assert "refstar" not in store.get_section("test_run")
    assert store.get_section("test_run")["teff"] == 4500
    request = setup.context()
    assert request is None or "refstar" not in request.config.preferences


def test_catalogue_edits_reach_the_dropdown(tmp_path):
    from refstarspage import RefStarsPage
    setup = _setup(tmp_path, _marcs_folder(tmp_path))
    page = RefStarsPage(on_change=setup.refresh_refstars)
    page.actions.value = "add  1"
    new = page.stars[-1]
    page.fields_input.value = f"{new.key}|label=Arcturus 2"
    page.fields_input.value = f"{new.key}|teff=4286 3"
    page.fields_input.value = f"{new.key}|logg=1.64 4"
    page.fields_input.value = f"{new.key}|feh=-0.52 5"
    page.fields_input.value = f"{new.key}|vmic=1.7 6"
    assert refstars.get(new.key).label == "Arcturus"
    assert dict(setup.run_form.fields["refstar"].options)[new.key] == "Arcturus"
    stored = store.get_section("reference_stars")["stars"]
    assert stored[-1] == {"key": new.key, "label": "Arcturus", "teff": 4286.0,
                          "logg": 1.64, "feh": -0.52, "vmic": 1.7}
    assert RefStarsPage().stars[-1].label == "Arcturus"     # a new session reads it back

    setup.run_form._on_input("value", "", f"refstar={new.key} 7")
    page.actions.value = f"del {new.key} 8"
    assert new.key not in refstars.by_key()
    assert setup.run_form.value("refstar") == refstars.CUSTOM_KEY   # the numbers stay
    assert setup.run_form.value("teff") == 4286

    page.actions.value = "reset  9"
    assert [s.key for s in refstars.stars()] == [s.key for s in refstars.DEFAULT_STARS]
    page.fields_input.value = "gdwarf|teff=1 10"
    assert refstars.DEFAULT_STARS[0].teff == 5750         # the defaults are copied, not shared


def test_grid_models_are_picked_from_a_filtered_list(tmp_path):
    folder = _marcs_folder(tmp_path)
    setup = _setup(tmp_path, folder)
    picker = setup.atmos_picker
    assert picker.node_mode
    assert "Pick model" in setup.grid_ui_html()
    labels = [label for label, *_ in picker._rows]
    assert labels == ["3800 K · log g 0.00 · [Fe/H] +0.00 · ξ 5",
                      "4500 K · log g 2.00 · [Fe/H] -1.00 · ξ 2",
                      "5750 K · log g 4.50 · [Fe/H] +0.00 · ξ 1",
                      "5750 K · log g 4.50 · [Fe/H] +0.00 · ξ 2"]
    picker.query.value_input = "5750 ξ 2"
    assert [label for label, *_ in picker._rows] == [labels[3]]
    picker.pick.value = f"{picker._rows[0][1]}\x00nonce"
    assert setup.stellar_values() == {"teff": 5750.0, "logg": 4.5, "feh": 0.0, "vmic": 2.0}
    assert picker.query.value == labels[3]
    assert setup.atmos_file() == str(folder / MARCS_NAMES[1])
    assert setup.run_form.value("refstar") == refstars.CUSTOM_KEY

    # the Individual atmosphere browses a folder again
    setup.grid_actions.value = "grid __individual__ 1"
    assert not picker.node_mode and picker.query.placeholder == "Type a path…"

    # an interpolated setup never scans the grid; entering pick mode does
    store.set_section("other_run", {"atmos_mode": "interp"})
    from synthsetup import SynthSetup
    lazy = SynthSetup("lazy", "other_run", setup.grids, grid_section="synthesis")
    assert lazy._nodes_stale and not lazy.atmos_picker._nodes
    lazy.grid_actions.value = "mode pick 1"
    assert not lazy._nodes_stale and len(lazy.atmos_picker._nodes) == 4


def test_strength_atmosphere_is_the_nearest_marcs_model(tmp_path):
    folder = _marcs_folder(tmp_path)
    strengths = linestrength.LineStrengths(master_path=tmp_path / "none",
                                           linelist=tmp_path / "none")
    assert not strengths.computable
    linestrength.set_marcs_grid(folder)
    assert strengths.computable
    assert linestrength.marcs_node(5700, 4.4, 0.1, 1.0).name == MARCS_NAMES[0]
    assert linestrength.marcs_node(5700, 4.4, 0.1, 1.8).name == MARCS_NAMES[1]
    assert linestrength.marcs_node(4400, 1.5, -1.0, 2.0).name == MARCS_NAMES[2]
    assert linestrength.marcs_node(*refstars.get("rsg").params.values()).name == MARCS_NAMES[3]
    linestrength.set_marcs_grid(tmp_path / "missing")
    assert not strengths.computable
    with pytest.raises(LookupError):
        linestrength.marcs_node(5750, 4.5, 0.0, 1.0)
