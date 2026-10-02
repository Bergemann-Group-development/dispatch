"""Namelist presets: lookup, stacking and where they sit among a run's settings.

What is load-bearing: a user file replaces the shipped preset of its name; later presets
win and the conflict is reported; the atmosphere selects the geometry preset; the request's
own settings, a grid's mesh and ``config.expert`` override presets; and a request preset
refines an NLTE step rather than being overridden by it.
"""
import logging

import pytest

from tumult.namelist import NameList
from tumult.synthesis import Snapshot, SynthConfig, SynthRequest, presets
from tumult.synthesis.workflows import preset_settings

STAR = (5777.0, 4.44, 0.0, 1.0)


@pytest.fixture(autouse=True)
def user_dir(tmp_path, monkeypatch):
    folder = tmp_path / "user_presets"
    folder.mkdir()
    monkeypatch.setattr(presets, "USER_DIR", folder)
    return folder


def test_user_preset_replaces_the_shipped_one(user_dir):
    (user_dir / "3d.nml").write_text("! My 3D mesh.\n&atmos_params nx=12 ny=12 nz=64 /\n")
    item = presets.load("3d")
    assert item.user and item.description == "My 3D mesh."
    assert item.settings == {"atmos_params": {"nx": 12, "ny": 12, "nz": 64}}
    assert [p.name for p in presets.available()].count("3d") == 1


def test_later_presets_win_and_the_conflict_is_reported(user_dir, caplog):
    (user_dir / "coarse.nml").write_text("&atmos_params nz=48 /\n")
    with caplog.at_level(logging.WARNING, logger=presets.__name__):
        settings, conflicts = presets.stack(["1d", "coarse"])
        presets.report(conflicts)
    assert settings["atmos_params"]["nz"] == 48
    assert settings["output_params"] == {"rt_to_file": False}
    assert conflicts == [("atmos_params", "nz", [("1d", 64), ("coarse", 48)])]
    assert "coarse wins" in caplog.text


def test_the_atmosphere_selects_the_geometry_preset(tmp_path):
    config = SynthConfig(root=tmp_path)
    one_d = SynthRequest(config=config, atmosphere="a.mod", star=STAR)
    assert preset_settings(one_d)["atmos_params"] == {"nz": 64}
    three_d = SynthRequest(config=config, atmosphere=Snapshot("cube", "M3DIS"), star=STAR)
    settings = preset_settings(three_d)
    assert settings["atmos_params"] == {"nx": 30, "ny": 30, "nz": 96}
    assert settings["run_params"] == {"rotate_atmos": True}
    assert all(settings[s] == {"sample_back_opac": True, "sample_linelist": True}
               for s in ("rt_iteration_params", "rt_output_params"))


def test_request_settings_and_expert_override_presets(tmp_path, user_dir):
    from tumult.synthesis.workflows import resolve_request
    (user_dir / "wide.nml").write_text(
        "&spectrum_params R_sampling=1000 /\n&linelist_params tolerance=1e-3 dlam=50 /\n")
    config = SynthConfig(root=tmp_path, resolving_power=50_000,
                         expert={"linelist_params": {"dlam": 20}})
    request = SynthRequest(config=config, atmosphere="a.mod", star=STAR,
                           windows=[(5000, 5010)], presets=["wide"])
    settings = resolve_request(request).settings
    assert settings["linelist_params"]["tolerance"] == 1e-3
    assert settings["linelist_params"]["dlam"] == 20
    assert settings["atmos_params"]["nz"] == 64


def test_a_request_preset_refines_the_departure_solve(tmp_path, user_dir, monkeypatch):
    from tumult.synthesis import api, workflows
    (user_dir / "loose.nml").write_text("&nlte_params convlim=1e-2 /\n")
    submitted = []

    class Run:
        has_spectrum, mode, folder = True, "NLTE", tmp_path

    node = (((5750.0, 4.5, 0.0), 1.0), 5750.0, 4.5, 0.0, 1.0, 1.0, "a.mod", "")
    monkeypatch.setattr(api, "select_nodes", lambda *a, **k: ([node], 1.0))
    monkeypatch.setattr(api, "_interp_crop", lambda *a: "spectrum")
    monkeypatch.setattr(api, "_run_batch", lambda jobs, *a, **k: (
        submitted.append(list(jobs)) or {key: Run() for key, _path, _nml in jobs}))
    request = SynthRequest(config=SynthConfig(root=tmp_path), atmosphere="a.mod",
                           star=STAR, presets=["loose"])
    workflows.compute(request, {"mode": "NLTE", "atoms": ["atom.mg"]}, set_spectrum=False)
    solve, spectrum = (nml for batch in submitted for _key, _path, nml in batch)
    assert solve["nlte_params"]["convlim"] == 1e-2
    assert solve["rt_output_params"]["angle_quad"] == "none"
    assert solve["atmos_params"]["save_atmos"] is True
    assert spectrum["atmos_params"]["read_departures"] is True
    assert spectrum["rt_iteration_params"]["maxiter"] == 0


def test_step_presets_are_chosen_by_the_departures(tmp_path):
    with pytest.raises(ValueError, match="follow from `departures`"):
        SynthRequest(config=SynthConfig(root=tmp_path), atmosphere="a.mod", star=STAR,
                     presets=["nlte"])


def test_a_comment_may_hold_an_apostrophe():
    text = "! the atom's own lines\n&linelist_params linelist_file='' dlam=5 /\n"
    assert NameList.from_string(text) == {"linelist_params": {"linelist_file": "",
                                                              "dlam": 5}}
