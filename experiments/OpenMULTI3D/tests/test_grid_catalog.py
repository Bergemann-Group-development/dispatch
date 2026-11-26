"""Model selection in a grid listing: exact parameters only, with the neighbours along one
parameter named when the node is missing."""
import pytest

from tumult.synthesis.grid_catalog import CATALOG, Listing, Model, RemoteFile


def _listing(*nodes):
    models = [Model(f"t{t:.0f}g{g * 10:.0f}m{abs(f) * 10:02.0f}", t, g, f,
                    {1: [RemoteFile(None, f"m{i}/atm3d", 100)]})
              for i, (t, g, f) in enumerate(nodes)]
    return Listing(CATALOG["stagger-dwarfs"], models)


def test_select_exact_node():
    listing = _listing((5000, 4.5, 0.0), (5777, 4.4, 0.0))
    assert listing.select(5777, 4.4, 0.0).teff == 5777


def test_select_names_neighbours_along_one_parameter():
    listing = _listing((5500, 4.5, 0.0), (6000, 4.5, 0.0), (6500, 4.5, 0.0),
                       (5777, 4.0, 0.0), (5777, 5.0, 0.0), (5777, 4.5, -1.0))
    with pytest.raises(LookupError) as err:
        listing.select(5777, 4.5, 0.0)
    text = str(err.value)
    assert "teff: t5500g45m00 (teff=5500), t6000g45m00 (teff=6000)" in text
    assert "logg: t5777g40m00 (logg=4), t5777g50m00 (logg=5)" in text
    assert "feh: t5777g45m10 (feh=-1)" in text
    assert "t6500" not in text


def test_select_without_any_neighbour():
    with pytest.raises(LookupError, match="none shares two"):
        _listing((5000, 4.5, 0.0)).select(6000, 4.0, -1.0)


def test_table_has_one_row_per_model():
    table = _listing((5000, 4.5, 0.0), (5777, 4.4, -0.5)).table()
    assert list(table.columns) == ["name", "teff", "logg", "feh", "snapshots", "size"]
    assert table["feh"].tolist() == [0.0, -0.5]
    assert table["size"].tolist() == [100, 100]
