"""Fetching published atmosphere grids in parts: the snapshot ordering that download and
selection share, the catalogue's reading of each archive, and the remote-zip transfer."""
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import threading
import zipfile

import pytest

from tumult.synthesis import grid_catalog, remote, snapshot_grid, workflows


# --- the snapshot ordering ------------------------------------------------------------------

def test_scrambled_halton_prefixes_are_stratified_and_nested():
    for m in range(1, 8):
        n = 2 ** m
        cells = sorted(int(snapshot_grid.scrambled_halton(i) * n) for i in range(n))
        assert cells == list(range(n))            # one point per dyadic interval
    for n in (1, 2, 3, 10, 11, 20):
        order = snapshot_grid.halton_order(n)
        assert sorted(order) == list(range(n))
        for k in range(1, n):
            assert set(order[:k]) < set(order[:k + 1])
    # a fixed seed: the same picks on every machine, so a download and a later
    # selection agree
    assert snapshot_grid.halton_order(10) == [7, 2, 9, 5, 4, 8, 0, 6, 3, 1]


def test_count_policy_selects_in_halton_order(tmp_path):
    node = tmp_path / "t50g45m00"
    node.mkdir()
    snaps = (8, 31, 32, 40, 51)
    for snap in snaps:
        (node / f"atm3d.t50g45m0005_{snap:05d}_80x80x240").write_text("")
    files = workflows.discover_atmos_grid(tmp_path, workflows.STAGGER_PATTERN, snapshots="2")
    picked = [int(p.name.split("_")[1]) for p in files[(5000., 4.5, 0.)]]
    order = snapshot_grid.halton_order(len(snaps))
    assert picked == sorted(snaps[i] for i in order[:2])
    assert picked == snapshot_grid.select_count(snaps, 2)


# --- the patterns of the catalogue grids ----------------------------------------------------

def test_stagger_pattern_reads_the_solar_models_too():
    p = workflows.compile_pattern(workflows.STAGGER_PATTERN)
    assert p.search("t5777g44m30/atm3d.t5777g44m3005_00008_80x80x240") == {
        "Teff": 5777., "logg": 4.4, "FeH": -3.0, "snap": 8.}
    assert p.search("t50g45m00/atm3d.t50g45m0005_00008_80x80x240")["Teff"] == 5000.


def test_m3dis_pattern_binds_the_folder_and_the_snapshot():
    p = workflows.compile_pattern(workflows.M3DIS_PATTERN)
    d = p.search("PLATO_t65.26g44.93m-3.000_v1.0/m3dis_999.atmos")
    assert (d["Teff"], d["FeH"], d["snap"]) == (6526., -3.0, 999.)
    assert abs(d["logg"] - 4.493) < 1e-9
    assert p.search("PLATO_t65.26g44.93m-3.000_v1.0/m3dis_999.mesh") is None


# --- a Stagger zip served locally -----------------------------------------------------------

def _stagger_zip(path):
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        for model, snaps in (("t50g45m00", (8, 31, 32)), ("t55g45m00", (3, 4, 5)),
                             ("t45g35m30", (1, 2)), ("t45g35m30_new", (7, 9))):
            for snap in snaps:
                zf.writestr(f"dwarfs/{model}/atm3d.{model.replace('_new', '')}05_{snap:05d}_80x80x240",
                            (f"cube {model} {snap} " * 500).encode())
                zf.writestr(f"dwarfs/{model}/mesh.{model.replace('_new', '')}05_{snap:05d}_80x80x240",
                            f"mesh {snap}".encode())
            zf.writestr(f"dwarfs/{model}/EOSrhoe.tab", b"eos")
        zf.writestr("dwarfs/README", b"readme")


class _RangeHandler(SimpleHTTPRequestHandler):
    """A file server honouring single byte ranges, as the grid hosts do."""

    def log_message(self, *args):
        pass

    def do_GET(self):
        header = self.headers.get("Range")
        path = Path(self.translate_path(self.path))
        if not header or not path.is_file():
            return super().do_GET()
        size = path.stat().st_size
        start, _, end = header.replace("bytes=", "").partition("-")
        start = int(start)
        end = int(end) if end else size - 1
        self.send_response(206)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.end_headers()
        with path.open("rb") as fh:
            fh.seek(start)
            self.wfile.write(fh.read(end - start + 1))


@pytest.fixture
def served(tmp_path):
    root = tmp_path / "www"
    root.mkdir()
    _stagger_zip(root / "dwarfs.zip")
    handler = partial(_RangeHandler, directory=str(root))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}/"
    server.shutdown()


def test_remote_zip_lists_and_extracts_members_by_range(served, tmp_path):
    archive = remote.RemoteZip(served + "dwarfs.zip")
    member = archive.members["dwarfs/t50g45m00/atm3d.t50g45m0005_00031_80x80x240"]
    assert member.method == 8 and member.usize == len("cube t50g45m00 31 ") * 500
    out = archive.extract(member, tmp_path / "out")
    assert out.read_bytes() == ("cube t50g45m00 31 " * 500).encode()
    assert archive.extract(member, out) == out          # complete copy kept


def test_plan_and_fetch_a_zenodo_grid_in_parts(served, tmp_path, monkeypatch):
    spec = grid_catalog.CATALOG["stagger-dwarfs"]
    local = grid_catalog.GridSpec(**{**spec.__dict__, "url": served + "dwarfs.zip"})
    monkeypatch.setitem(grid_catalog.CATALOG, "stagger-dwarfs", local)
    dest = tmp_path / "grid"
    listing = grid_catalog.list_models("stagger-dwarfs")
    # the superseded model is not listed; every listed model has its snapshots
    assert [m.name for m in listing.models] == ["t45g35m30_new", "t50g45m00", "t55g45m00"]
    assert listing.find("t50g45m00").available == [8, 31, 32]
    plan = grid_catalog.plan(listing, dest, snapshots=1, models=["t50g45m00", "t55g45m00"])
    assert all(len(v) == 1 for v in plan.chosen.values())
    assert len(plan.files) == 4 and not plan.present
    record = grid_catalog.fetch(plan)
    assert record == {"id": "stagger-dwarfs", "name": local.name, "folder": str(dest),
                      "pattern": workflows.STAGGER_PATTERN, "format": "Stagger2",
                      "snapshots": "1"}
    manifest = grid_catalog.read_manifest(dest)
    assert manifest["snapshots"] == 1 and set(manifest["fetched"]) == {"t50g45m00", "t55g45m00"}
    assert listing.find("t50g45m00").on_disk(dest) == manifest["fetched"]["t50g45m00"]
    # the fetched layout is what the registered record scans, and the run-time count
    # policy selects exactly the fetched snapshots
    files = snapshot_grid.discover_grid(record)
    assert {k: [p.name for p in v] for k, v in files.items()} == {
        (5000., 4.5, 0.): [f"atm3d.t50g45m0005_{manifest['fetched']['t50g45m00'][0]:05d}_80x80x240"],
        (5500., 4.5, 0.): [f"atm3d.t55g45m0005_{manifest['fetched']['t55g45m00'][0]:05d}_80x80x240"]}
    # more snapshots later only add files
    plan2 = grid_catalog.plan(listing, dest, snapshots=2)
    assert len(plan2.present) == 4 and len(plan2.pending) == 8
    grid_catalog.fetch(plan2)
    assert grid_catalog.read_manifest(dest)["snapshots"] == 2
    assert len(list(dest.glob("*/atm3d.*"))) == 6
    assert len(listing.find("t50g45m00").on_disk(dest)) == 2
