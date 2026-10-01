"""A local stand-in for the Seafile share that publishes the master line lists.

:func:`serve_share` runs a threaded HTTP server over a folder of files and answers the
two requests :mod:`tumult.synthesis.linelist_catalog` (and :mod:`.input_data`) makes of a share link: the
``dirents`` listing of a folder as Seafile's share-link API shapes it (``dirent_list`` of
``file_name`` / ``size`` / ``is_dir`` entries), and a file download through the
``/d/<token>/files/?p=<path>&dl=1`` form, served with single byte ranges as Keeper does.
"""
from __future__ import annotations

from contextlib import contextmanager
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import threading
import urllib.parse

TOKEN = "testshare"


class _ShareHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def _send_bytes(self, body, content_type):
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        url = urllib.parse.urlsplit(self.path)
        query = urllib.parse.parse_qs(url.query)
        root = Path(self.directory)
        if url.path == f"/api/v2.1/share-links/{TOKEN}/dirents/":
            folder = root / query["path"][0].strip("/")
            entries = []
            for item in sorted(folder.iterdir()):
                if item.is_dir():
                    entries.append({"folder_name": item.name, "is_dir": True,
                                    "folder_path": f"/{item.relative_to(root)}/"})
                else:
                    entries.append({"file_name": item.name, "is_dir": False, "size":
                                    item.stat().st_size, "file_path":
                                    f"/{item.relative_to(root)}"})
            self._send_bytes(json.dumps({"dirent_list": entries}).encode(),
                             "application/json")
            return
        if url.path == f"/d/{TOKEN}/files/":
            path = root / query["p"][0].strip("/")
        else:
            self.send_error(404)
            return
        if not path.is_file():
            self.send_error(404)
            return
        size = path.stat().st_size
        header = self.headers.get("Range")
        start, end = 0, size - 1
        if header:
            first, _, last = header.replace("bytes=", "").partition("-")
            start, end = int(first), int(last) if last else size - 1
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        else:
            self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(end - start + 1))
        self.end_headers()
        with path.open("rb") as fh:
            fh.seek(start)
            self.wfile.write(fh.read(end - start + 1))


@contextmanager
def serve_share(root):
    """Serve ``root`` as a share link; yields the ``ShareFolder`` fields ``server`` and
    ``token`` as a tuple ``(server_url, token)``."""
    handler = partial(_ShareHandler, directory=str(root))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_address[1]}", TOKEN
    finally:
        server.shutdown()


def write_lists(root, folder="linelists/master_linelists"):
    """Three small line-list files (and a sub-folder, which listings skip) under
    ``root/folder``; returns ``{name: content}``."""
    target = Path(root) / folder
    target.mkdir(parents=True)
    files = {"12C16O.txt": b"'C O' 3\n" + b"4300.0 1.0 -2.0\n" * 50,
             "Hlinedata": b"H lines\nnlines 1\n6562.8 2 3 10.2 12.1 0.71 'Ha' 4.7e8\n",
             "vald_atoms.txt": b"'Fe I' 100\n" + b"5000.0 3.0 -1.0\n" * 400}
    for name, body in files.items():
        (target / name).write_bytes(body)
    (target / "extra").mkdir()
    (target / "extra" / "README").write_bytes(b"not a line list")
    return files
