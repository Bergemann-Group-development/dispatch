"""The published line lists OpenMULTI3D can fetch, and how to fetch them.

The master line lists (one file per species, plus HOP's ``Hlinedata``) are hosted on the
MPCDF Keeper share that also carries the other large inputs (``README_INSTALL.md``);
:data:`SOURCE` names its ``linelists/master_linelists`` folder. A download is planned
before any bytes move: :func:`list_files` reads the share's folder once through the
Seafile share-link API, :func:`plan` takes the chosen files, checks what the destination
already holds and totals the transfer, and :func:`fetch` carries the plan out, one
resumable transfer per file. Downloading into the folder the synthesis reads its line
list from makes each list available as soon as it is complete: the readers take a
folder's direct children, and a transfer is a ``.part`` file until it is whole.

The command line is :mod:`.input_data`, which serves this and the other share folders.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import urllib.parse

from .grid_catalog import human
from . import remote


@dataclass(frozen=True)
class ShareFolder:
    """A folder of a public Seafile share link: ``server`` and share ``token`` address
    the share, ``path`` the folder within it. A ``nested`` folder is read with its
    sub-folders, and its files are named by their path below ``path``."""
    id: str
    name: str
    server: str
    token: str
    path: str
    citation: str
    nested: bool = False

    @property
    def page(self):
        """Where the share is browsed."""
        return f"{self.server}/d/{self.token}/?p={urllib.parse.quote(self.path.rstrip('/'))}"

    def listing_url(self, sub=""):
        """The listing of ``path``, or of its sub-folder ``sub`` (ending in ``/``)."""
        return (f"{self.server}/api/v2.1/share-links/{self.token}/dirents/?"
                + urllib.parse.urlencode({"path": self.path + sub},
                                         quote_via=urllib.parse.quote))

    def file_url(self, name):
        return (f"{self.server}/d/{self.token}/files/?"
                + urllib.parse.urlencode({"p": f"{self.path.rstrip('/')}/{name}", "dl": "1"},
                                         quote_via=urllib.parse.quote))


SOURCE = ShareFolder(
    id="keeper-master", name="Master line lists",
    server="https://keeper.mpdl.mpg.de", token="3191475aebe340638422",
    path="/linelists/master_linelists/",
    citation="OpenMULTI3D input data on MPCDF Keeper")


@dataclass(frozen=True)
class RemoteList:
    """One file of the share: its ``name`` under the folder (a relative path in a
    ``nested`` source), ``size`` in bytes as the share reports it, and the ``url`` it
    downloads from."""
    name: str
    size: int
    url: str

    def on_disk(self, dest):
        """Whether ``dest`` holds the complete file (a ``.part`` in progress does not)."""
        path = Path(dest) / self.name
        return path.is_file() and path.stat().st_size == self.size


@dataclass
class Listing:
    """The share folder as listed: its files sorted by name."""
    source: ShareFolder
    files: list

    def find(self, name):
        return next((f for f in self.files if f.name == name), None)

    @property
    def size(self):
        return sum(f.size for f in self.files)


def list_files(source=SOURCE, cancel=None):
    """Read the share folder once: its files with their sizes. Sub-folders are descended
    only for a ``nested`` source; the line-list readers take a folder's direct children.
    ``cancel`` is a :class:`~.cancel.CancelToken` checked before each request."""
    files = []
    folders = [""]
    while folders:
        sub = folders.pop()
        if cancel is not None:
            cancel.check()
        for e in remote.read_json(source.listing_url(sub))["dirent_list"]:
            if not e["is_dir"]:
                name = sub + e["file_name"]
                files.append(RemoteList(name, int(e["size"]), source.file_url(name)))
            elif source.nested:
                folders.append(f"{sub}{e['folder_name']}/")
    return Listing(source, sorted(files, key=lambda f: f.name))


@dataclass
class Plan:
    listing: Listing
    dest: Path
    files: list            # the chosen RemoteList, in transfer order
    present: list          # those already complete under dest

    @property
    def total_bytes(self):
        return sum(f.size for f in self.files)

    @property
    def present_bytes(self):
        return sum(f.size for f in self.present)

    @property
    def pending(self):
        done = {f.name for f in self.present}
        return [f for f in self.files if f.name not in done]

    def summary(self):
        return "\n".join([
            f"{self.listing.source.name}: {self.listing.source.page}",
            f"destination: {self.dest}",
            f"{human(self.total_bytes)} in {len(self.files)} files, "
            f"{human(self.present_bytes)} already present"])


def plan(listing, dest, names=None):
    """Plan fetching ``names`` (all of the listing by default) into ``dest``."""
    dest = Path(dest).expanduser()
    files = listing.files if names is None else [
        f for f in listing.files if f.name in set(names)]
    return Plan(listing, dest, files, [f for f in files if f.on_disk(dest)])


def fetch(plan, progress=None, cancel=None):
    """Carry out ``plan``: fetch every pending file. ``progress(done_bytes, total_bytes,
    name)`` is called as bytes arrive; ``cancel`` is a :class:`~.cancel.CancelToken`."""
    plan.dest.mkdir(parents=True, exist_ok=True)
    total = plan.total_bytes
    done = plan.present_bytes
    for f in plan.pending:
        if cancel is not None:
            cancel.check()

        def report(got, _size, base=done, name=f.name):
            if progress is not None:
                progress(base + got, total, name)
        remote.fetch_file(f.url, plan.dest / f.name, size=f.size, progress=report,
                          cancel=cancel)
        done += f.size
    return plan.dest


def delete(listing, dest, names=None):
    """Remove the listed files ``names`` (all of the listing by default) from ``dest``,
    with any partial transfers of theirs. Files the share does not list are the user's
    own and stay."""
    dest = Path(dest).expanduser()
    for f in listing.files if names is None else [f for f in listing.files
                                                    if f.name in set(names)]:
        (dest / f.name).unlink(missing_ok=True)
        (dest / (f.name + ".part")).unlink(missing_ok=True)
