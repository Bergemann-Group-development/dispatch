"""The line-list panel marks exactly the picked transition, identified by wavelength and ion."""
from types import SimpleNamespace

import pandas as pd
import pytest

from linenav import LAM0_UNIT, LineNavigator


@pytest.fixture
def nav():
    # the row-marking logic needs only the focus and the display unit, not a document
    nav = object.__new__(LineNavigator)
    nav.canvas = SimpleNamespace(xunit=LAM0_UNIT)
    nav._scroll = "none"
    return nav


@pytest.fixture
def rows():
    # a Fe I line with a TiO line 0.4 mÅ away and a second Fe I component 0.6 mÅ away
    return pd.DataFrame({
        "ion": ["TiO", "Fe I", "Fe I", "Fe II"],
        "lam0": [5000.0996, 5000.1000, 5000.1006, 5010.0],
        "evi": [0.1, 2.0, 2.0, 3.0], "loggf": [-1.0, -0.5, -2.0, -1.0],
        "strength": [1.0, 5.0, 0.1, 2.0],
    })


def marked(nav, rows):
    return nav._body_html(rows).count("jump-target")


def test_a_pick_marks_only_its_own_ion_and_the_nearest_line(nav, rows):
    nav.focus = (5000.1000, "Fe I")
    assert nav._focus_row(rows) == 1
    assert marked(nav, rows) == 1
    nav.focus = (5000.0996, "TiO")
    assert nav._focus_row(rows) == 0


def test_a_row_click_carries_its_ion_through_the_action_value(nav, rows, monkeypatch):
    main = SimpleNamespace(x_range=SimpleNamespace(start=0.0, end=1.0, bounds=None))
    nav.canvas = SimpleNamespace(xunit=LAM0_UNIT, main=main,
                                 nav=SimpleNamespace(x_range=SimpleNamespace(
                                     start=0.0, end=1.0, bounds=None)))
    nav.jump_distance = 10.0
    monkeypatch.setattr(nav, "_flash", lambda x: None, raising=False)
    # the page sends the ion URI-encoded, as encodeURIComponent("Fe I") does
    nav.jump("5000.1000|Fe%20I")
    assert nav.focus == (5000.1, "Fe I")
    assert nav._focus_row(rows) == 1


def test_no_row_is_marked_when_the_picked_ion_is_absent(nav, rows):
    nav.focus = (5000.1000, "Cr I")
    assert nav._focus_row(rows) is None
    assert marked(nav, rows) == 0
