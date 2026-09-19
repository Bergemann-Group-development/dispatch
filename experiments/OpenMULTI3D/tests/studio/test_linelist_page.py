"""The Line lists preferences page.

What matters: the page lists the share once when opened, shows per file whether it is
on disk, downloads exactly the file (or all files) the user asks for into the folder
the *Line list* preference names, reports the outcome, and tells the document that
folder changed so the master list is re-read. When the preference names a single file
the page says where downloads go instead. The share is the local stand-in of
:mod:`tests.support.share`; the IO pool is replaced by a context that runs each job
at once, so the list / fetch / settle sequence is followed inline.
"""
import re

from tests.support.share import serve_share, write_lists

_FETCH_BTN = re.compile(r'linelist-fetch" data-name="([^"]+)"')
_DEL_BTN = re.compile(r'linelist-del" data-name="([^"]+)"')


class _SyncContext:
    def __init__(self):
        self.notes = []
        self.jobs = []

    def notify(self, message, level="warn"):
        self.notes.append(message)

    def submit_io(self, name, work, done):
        from tumult.synthesis.cancel import CancelToken
        self.jobs.append(name)
        done(work(CancelToken()))
        return name

    def schedule(self, fn):
        fn()


def _manager(server, token, ctx, changes):
    from tumult.synthesis import linelist_catalog
    from TumultStudio.preferences.linelistmanager import LineListManager
    source = linelist_catalog.ShareFolder(
        id="test", name="Test lists", server=server, token=token,
        path="/linelists/master_linelists/", citation="test")
    return LineListManager(on_change=changes.append, ctx=ctx, source=source)


def test_page_lists_downloads_and_deletes_into_the_preference_folder(tmp_path, studio_config):
    from TumultStudio.workspace import store
    files = write_lists(tmp_path / "www")
    dest = tmp_path / "work" / "lists"
    store.set_section("synthesis", {"dispatch_cwd": str(tmp_path / "work"),
                                    "linelist": "lists"})
    ctx, changes = _SyncContext(), []
    with serve_share(tmp_path / "www") as (server, token):
        lm = _manager(server, token, ctx, changes)
        page = lm.list_html()
        assert f"Downloads go to <b>{dest}</b>, the line list of the Synthesis page" in page
        assert "not listed" in page and "linelist-table" not in page

        # opening the page lists the share once; every file offers a download
        lm.actions.value = "list 1"
        lm.actions.value = "list 2"
        assert ctx.jobs == ["list-lines Test lists"]
        page = lm.list_html()
        assert "0 of 3 files on disk" in page and "linelist-fetch-all" in page
        assert _FETCH_BTN.findall(page) == sorted(files) and not _DEL_BTN.findall(page)
        assert "linelist-del-all" not in page

        # one file: only it arrives, its row flips to delete, the folder is reported
        lm.actions.value = "fetch Hlinedata 3"
        assert ctx.notes == ["Hlinedata downloaded"]
        assert [p.name for p in dest.iterdir()] == ["Hlinedata"]
        assert (dest / "Hlinedata").read_bytes() == files["Hlinedata"]
        assert changes == [dest]
        page = lm.list_html()
        assert "1 of 3 files on disk" in page and "grid-progress" not in page
        assert _DEL_BTN.findall(page) == ["Hlinedata"]
        assert _FETCH_BTN.findall(page) == ["12C16O.txt", "vald_atoms.txt"]
        lm.actions.value = "fetch Hlinedata 4"
        assert ctx.notes[-1] == "Hlinedata already on disk"

        # sorting by size puts the largest first
        lm.actions.value = "sort size 5"
        names = re.findall(r'class="line-ion" title="([^"]+)"', lm.list_html())
        assert names == ["vald_atoms.txt", "12C16O.txt", "Hlinedata"]

        # all files, then delete one and all: the user's own file stays
        lm.actions.value = "fetch-all 6"
        assert ctx.notes[-1] == "test lists downloaded"
        assert sorted(p.name for p in dest.iterdir()) == sorted(files)
        page = lm.list_html()
        assert "3 of 3 files on disk" in page and "linelist-fetch-all" not in page
        (dest / "mine.txt").write_bytes(b"'Fe I' 1\n")
        lm.actions.value = "del 12C16O.txt 7"
        assert not (dest / "12C16O.txt").exists() and changes[-1] == dest
        lm.actions.value = "del-all 8"
        assert [p.name for p in dest.iterdir()] == ["mine.txt"]
        assert "0 of 3 files on disk" in lm.list_html()


def test_page_names_the_fallback_folder_when_the_preference_is_a_file(tmp_path, studio_config):
    from TumultStudio.workspace import store
    from TumultStudio.preferences.linelistmanager import LineListManager
    (tmp_path / "sun.txt").write_text("'Fe I' 1\n")
    store.set_section("synthesis", {"dispatch_cwd": str(tmp_path), "linelist": "sun.txt"})
    lm = LineListManager()
    page = lm.list_html()
    fallback = tmp_path / "input" / "linelists" / "master_linelists"
    assert f"Downloads go to <b>{fallback}</b>. The line list of the Synthesis page names a single file" in page
    assert lm.folder == fallback


def test_listing_failure_is_shown_with_a_retry(tmp_path, studio_config, monkeypatch):
    from tumult.synthesis import remote

    class _Ctx(_SyncContext):
        def submit_io(self, name, work, done):
            from tumult.synthesis.cancel import CancelToken
            import pytest
            with pytest.raises(remote.Unreachable):
                work(CancelToken())
            return name
    monkeypatch.setattr(remote, "LIST_TIMEOUT", 1)
    monkeypatch.setattr(remote, "LIST_RETRIES", 1)
    ctx = _Ctx()
    lm = _manager("http://127.0.0.1:9", "dead", ctx, [])
    lm.actions.value = "list 1"
    page = lm.list_html()
    assert "could not list the share: 127.0.0.1 is not responding" in page
    assert "linelist-relist" in page and ctx.notes[0].startswith("Could not list Test lists")
