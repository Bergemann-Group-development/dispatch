"""tumult.downloads: one interface to the grids and the datasets, on made-up listings so
that nothing is downloaded."""
import pytest

from tumult import downloads
from tumult.synthesis import grid_catalog, input_data, linelist_catalog
from tumult.synthesis.grid_catalog import Model, RemoteFile
from tumult.synthesis.linelist_catalog import RemoteList
from tests.support.share import serve_share, write_lists


@pytest.fixture
def sources(monkeypatch):
    """A one-model Stagger listing and a two-file atoms listing; transfers write the files."""
    model = Model("t5777g44m00", 5777, 4.4, 0.0, {
        snap: [RemoteFile(None, f"t5777g44m00/{kind}.t5777g44m0005_{snap:05d}_80x80x240", 3)
               for kind in ("atm3d", "mesh")] for snap in (76, 459)})
    listings = {
        "stagger-dwarfs": grid_catalog.Listing(grid_catalog.CATALOG["stagger-dwarfs"], [model]),
        "atoms": linelist_catalog.Listing(input_data.DATASETS["atoms"].source,
                                          [RemoteList("atom.ba06", 3, ""),
                                           RemoteList("atom.mg86g", 3, "")]),
    }
    monkeypatch.setattr(downloads, "_listing", listings.__getitem__)
    transfers = []

    def write(todo, progress=None):
        transfers.append(todo)
        for f in todo.pending:
            path = todo.dest / getattr(f, "path", getattr(f, "name", None))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"abc")
    monkeypatch.setattr(grid_catalog, "fetch", write)
    monkeypatch.setattr(linelist_catalog, "fetch", write)
    return transfers


def test_grid_model_by_parameters(sources, tmp_path):
    paths = downloads.fetch("stagger-dwarfs", teff=5777, logg=4.4, feh=0.0, dest=tmp_path)
    assert [p.name for p in paths] == ["atm3d.t5777g44m0005_00459_80x80x240"]
    assert all(p.is_file() for p in paths)
    assert downloads.table("stagger-dwarfs", dest=tmp_path)["on_disk"].tolist() == [1]


def test_dataset_files_by_name(sources, tmp_path):
    paths = downloads.fetch("atoms", "atom.mg86g", dest=tmp_path)
    assert paths == [tmp_path / "atom.mg86g"]
    assert downloads.table("atoms", dest=tmp_path)["on_disk"].tolist() == [False, True]


def test_complete_files_are_not_fetched_again(sources, tmp_path, capsys):
    downloads.fetch("atoms", "atom.ba06", dest=tmp_path)
    downloads.fetch("atoms", "atom.ba06", dest=tmp_path, verbose=True)
    assert len(sources) == 1
    assert "already on disk" in capsys.readouterr().out


@pytest.mark.parametrize("call, error", [
    (lambda: downloads.fetch("atoms", teff=5000), TypeError),
    (lambda: downloads.fetch("stagger-dwarfs", "atom.ba06"), TypeError),
    (lambda: downloads.fetch("stagger-dwarfs", teff=5777, logg=4.4), TypeError),
    (lambda: downloads.fetch("atoms", "atom.xx"), LookupError),
    (lambda: downloads.fetch("stagger-dwarfs", teff=5800, logg=4.4, feh=0.0), LookupError),
    (lambda: downloads.fetch("nonexistent"), LookupError),
])
def test_wrong_selection_raises_before_any_transfer(sources, call, error):
    with pytest.raises(error):
        call()
    assert sources == []


def test_command_line_fetches_a_dataset_into_the_chosen_folder(tmp_path, capsys, monkeypatch):
    files = write_lists(tmp_path / "www")
    dest = tmp_path / "dest"
    with serve_share(tmp_path / "www") as (server, token):
        source = linelist_catalog.ShareFolder(
            id="test", name="Test lists", server=server, token=token,
            path="/linelists/master_linelists/", citation="test")
        monkeypatch.setitem(input_data.DATASETS, "atoms",
                            input_data.Dataset("atoms", source, "atoms/public"))
        downloads._listing.cache_clear()
        assert downloads.main(["list", "atoms"]) == 0
        assert "Hlinedata" in capsys.readouterr().out
        assert downloads.main(["fetch", "atoms", "--dest", str(dest), "--file", "Hlinedata"]) == 0
        assert [p.name for p in dest.iterdir()] == ["Hlinedata"]
        assert (dest / "Hlinedata").read_bytes() == files["Hlinedata"]
    downloads._listing.cache_clear()
