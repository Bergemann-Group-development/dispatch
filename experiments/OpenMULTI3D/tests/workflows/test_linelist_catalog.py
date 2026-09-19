"""Fetching the published master line lists.

What matters: a user chooses individual line lists from the share and gets exactly
those files, complete and unchanged, into the folder the synthesis reads; a
re-run moves no bytes for files already there; deleting touches nothing the share
does not list. The share is stood in for by :mod:`tests.support.share`, which answers
the same two requests as Keeper (the dirents listing and a ranged file download), so a
mismatch here means the client side of that protocol changed.
"""
from tumult.synthesis import linelist_catalog
from tests.support.share import serve_share, write_lists


def _source(server, token):
    return linelist_catalog.ShareFolder(
        id="test", name="Test lists", server=server, token=token,
        path="/linelists/master_linelists/", citation="test")


def test_listing_plan_fetch_and_delete_follow_the_chosen_files(tmp_path):
    files = write_lists(tmp_path / "www")
    dest = tmp_path / "dest"
    with serve_share(tmp_path / "www") as (server, token):
        source = _source(server, token)
        listing = linelist_catalog.list_files(source)
        # files only, by name, with the sizes the share reports
        assert [f.name for f in listing.files] == sorted(files)
        assert listing.find("Hlinedata").size == len(files["Hlinedata"])
        assert listing.size == sum(len(b) for b in files.values())

        chosen = linelist_catalog.plan(listing, dest, ["Hlinedata", "12C16O.txt"])
        assert [f.name for f in chosen.pending] == ["12C16O.txt", "Hlinedata"]
        assert chosen.total_bytes == len(files["Hlinedata"]) + len(files["12C16O.txt"])
        seen = []
        linelist_catalog.fetch(chosen, progress=lambda done, total, name: seen.append(
            (done, total, name)))
        assert sorted(p.name for p in dest.iterdir()) == ["12C16O.txt", "Hlinedata"]
        assert (dest / "Hlinedata").read_bytes() == files["Hlinedata"]
        # progress counts bytes cumulatively over the plan and ends at the total
        assert seen[-1][:2] == (chosen.total_bytes, chosen.total_bytes)
        assert all(a[0] <= b[0] for a, b in zip(seen, seen[1:]))

        # the whole listing later: the two present files are skipped, one is fetched
        again = linelist_catalog.plan(listing, dest)
        assert [f.name for f in again.present] == ["12C16O.txt", "Hlinedata"]
        assert [f.name for f in again.pending] == ["vald_atoms.txt"]
        linelist_catalog.fetch(again)
        assert all(f.on_disk(dest) for f in listing.files)

        # a partial transfer is not "on disk"; deleting takes it along, and leaves the
        # user's own file alone
        (dest / "Hlinedata.part").write_bytes(b"x")
        (dest / "Hlinedata").unlink()
        (dest / "mine.txt").write_bytes(b"'Fe I' 1\n")
        assert not listing.find("Hlinedata").on_disk(dest)
        linelist_catalog.delete(listing, dest, ["Hlinedata", "vald_atoms.txt"])
        assert sorted(p.name for p in dest.iterdir()) == ["12C16O.txt", "mine.txt"]
        linelist_catalog.delete(listing, dest)
        assert [p.name for p in dest.iterdir()] == ["mine.txt"]


def test_share_urls_address_the_folder_and_its_files():
    source = linelist_catalog.SOURCE
    assert source.listing_url() == ("https://keeper.mpdl.mpg.de/api/v2.1/share-links/"
                                    "3191475aebe340638422/dirents/"
                                    "?path=%2Flinelists%2Fmaster_linelists%2F")
    assert source.file_url("Hlinedata") == (
        "https://keeper.mpdl.mpg.de/d/3191475aebe340638422/files/"
        "?p=%2Flinelists%2Fmaster_linelists%2FHlinedata&dl=1")
    assert source.page.startswith("https://keeper.mpdl.mpg.de/d/3191475aebe340638422/")
