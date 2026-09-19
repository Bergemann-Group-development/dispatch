"""HTTP access to the archives that publish model-atmosphere grids and line lists.

The grids in :mod:`.grid_catalog` come as one large zip per grid (Zenodo), one tar of
gzipped models (the MARCS server) or a web directory index (the M3DIS server); the line
lists of :mod:`.linelist_catalog` are plain files on a Seafile share listed through its
JSON API. Every one of those hosts honours byte-range requests, so an archive can be
fetched in parts: a zip's central directory names each member's offset and a member is
inflated straight from its range, an index is read from the HTML listing, and a plain
file resumes from its partial copy. Only the standard library is used.

``progress(done, total)`` callbacks count bytes of the current transfer; ``cancel`` is a
:class:`~tumult.synthesis.cancel.CancelToken` checked between chunks. Listing a grid
(directory pages, a zip's central directory, a HEAD) uses the short ``LIST_TIMEOUT`` with
one retry, so an unreachable host fails within a minute; transfers keep the patient
``TIMEOUT`` and ``RETRIES``.
"""
from __future__ import annotations

from dataclasses import dataclass
import html
import json
from pathlib import Path
import re
import ssl
import struct
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib

CHUNK = 1 << 20
RETRIES = 3
TIMEOUT = 120
LIST_TIMEOUT = 20
LIST_RETRIES = 2


def _ssl_context():
    """The python.org macOS builds ship without a CA trust store; the certifi bundle
    (installed with the project's dependencies) fills that gap when present."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


_SSL = _ssl_context()


class Unreachable(RuntimeError):
    """The host did not answer a request within its attempts."""


def _open(url, headers=None, method=None, quick=False):
    """An open response for ``url``; ``quick`` uses the listing timeout and retries."""
    req = urllib.request.Request(url, headers=headers or {}, method=method)
    timeout, retries = (LIST_TIMEOUT, LIST_RETRIES) if quick else (TIMEOUT, RETRIES)
    last = None
    for attempt in range(1, retries + 1):
        try:
            return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
        except urllib.error.HTTPError as exc:
            if exc.code < 500:
                raise
            last = exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
        if attempt < retries:
            time.sleep(2 * attempt)
    host = urllib.parse.urlsplit(url).hostname or url
    raise Unreachable(f"{host} is not responding ({last})")


def read_json(url):
    """The JSON document at ``url``, with the listing timeout."""
    with _open(url, headers={"Accept": "application/json"}, quick=True) as resp:
        return json.loads(resp.read())


def content_length(url):
    """The size in bytes the server reports for ``url``."""
    with _open(url, method="HEAD", quick=True) as resp:
        return int(resp.headers["Content-Length"])


def fetch_range(url, start, end):
    """Bytes ``start..end`` (inclusive) of ``url``, with the listing timeout: the callers
    read archive directories and headers, never bulk data."""
    with _open(url, headers={"Range": f"bytes={start}-{end}"}, quick=True) as resp:
        if resp.status != 206:
            raise RuntimeError(f"{url} does not serve byte ranges")
        return resp.read()


def _copy(resp, fh, done, total, progress, cancel):
    while True:
        if cancel is not None:
            cancel.check()
        buf = resp.read(CHUNK)
        if not buf:
            return done
        fh.write(buf)
        done += len(buf)
        if progress is not None:
            progress(done, total)


def fetch_file(url, dest, size=None, progress=None, cancel=None):
    """Download ``url`` to ``dest``, resuming a ``.part`` left by an earlier attempt.

    ``size`` (the expected length) lets an existing complete file be skipped without a
    request; it is read from the server otherwise."""
    dest = Path(dest)
    if size is None:
        size = content_length(url)
    if dest.exists() and dest.stat().st_size == size:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    if have >= size:
        have = 0
    headers = {"Range": f"bytes={have}-"} if have else {}
    with _open(url, headers=headers) as resp:
        if have and resp.status != 206:
            have = 0
        with part.open("ab" if have else "wb") as fh:
            done = _copy(resp, fh, have, size, progress, cancel)
    if done != size:
        raise RuntimeError(f"{url}: got {done} of {size} bytes")
    part.replace(dest)
    return dest


# --- zip archives read in place ------------------------------------------------------

@dataclass(frozen=True)
class ZipMember:
    name: str
    offset: int      # of the local file header
    csize: int
    usize: int
    method: int      # 0 stored, 8 deflate
    crc: int


def _zip64(extra, *fields):
    """Replace the 0xFFFFFFFF ``fields`` from the zip64 extra record, in its order."""
    values = list(fields)
    q = 0
    while q + 4 <= len(extra):
        hid, hsz = struct.unpack("<HH", extra[q:q + 4])
        if hid == 1:
            wide = list(struct.unpack("<" + "Q" * (hsz // 8), extra[q + 4:q + 4 + hsz]))
            values = [wide.pop(0) if v == 0xFFFFFFFF and wide else v for v in values]
            break
        q += 4 + hsz
    return values


def parse_central_directory(data):
    """The members listed in a zip central directory ``data`` (zip64 aware)."""
    members = []
    p = 0
    while p + 46 <= len(data) and data[p:p + 4] == b"PK\x01\x02":
        method, = struct.unpack("<H", data[p + 10:p + 12])
        crc, csize, usize, nlen, xlen, clen = struct.unpack("<IIIHHH", data[p + 16:p + 34])
        offset, = struct.unpack("<I", data[p + 42:p + 46])
        name = data[p + 46:p + 46 + nlen].decode("utf-8")
        extra = data[p + 46 + nlen:p + 46 + nlen + xlen]
        usize, csize, offset = _zip64(extra, usize, csize, offset)
        members.append(ZipMember(name, offset, csize, usize, method, crc))
        p += 46 + nlen + xlen + clen
    return members


class RemoteZip:
    """A zip served over HTTP, listed from its central directory and read per member.

    Reading the directory costs two or three small range requests; each member is then
    one request for its compressed bytes, inflated on the way to disk."""

    def __init__(self, url):
        self.url = url
        self.size = content_length(url)
        self.members = {m.name: m for m in parse_central_directory(self._directory())}

    def _directory(self):
        tail = fetch_range(self.url, max(0, self.size - 65536), self.size - 1)
        i = tail.rfind(b"PK\x05\x06")
        if i < 0:
            raise RuntimeError(f"{self.url} is not a zip archive")
        cd_size, cd_off = struct.unpack("<II", tail[i + 12:i + 20])
        if 0xFFFFFFFF in (cd_size, cd_off):
            j = tail.rfind(b"PK\x06\x07")
            z64_off, = struct.unpack("<Q", tail[j + 8:j + 16])
            z64 = fetch_range(self.url, z64_off, z64_off + 55)
            cd_size, cd_off = struct.unpack("<QQ", z64[40:56])
        return fetch_range(self.url, cd_off, cd_off + cd_size - 1)

    def extract(self, member, dest, progress=None, cancel=None):
        """Write ``member`` to ``dest``; a complete copy of the right size is kept."""
        dest = Path(dest)
        if dest.exists() and dest.stat().st_size == member.usize:
            return dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        head = fetch_range(self.url, member.offset, member.offset + 29)
        if head[:4] != b"PK\x03\x04":
            raise RuntimeError(f"{self.url}: bad local header for {member.name}")
        nlen, xlen = struct.unpack("<HH", head[26:30])
        start = member.offset + 30 + nlen + xlen
        inflate = zlib.decompressobj(-15) if member.method == 8 else None
        if member.method not in (0, 8):
            raise RuntimeError(f"{member.name}: unsupported compression {member.method}")
        part = dest.with_name(dest.name + ".part")
        done = 0
        crc = 0
        with _open(self.url, headers={"Range": f"bytes={start}-{start + member.csize - 1}"}) \
                as resp, part.open("wb") as fh:
            while True:
                if cancel is not None:
                    cancel.check()
                buf = resp.read(CHUNK)
                if not buf:
                    break
                done += len(buf)
                out = inflate.decompress(buf) if inflate else buf
                crc = zlib.crc32(out, crc)
                fh.write(out)
                if progress is not None:
                    progress(done, member.csize)
            if inflate:
                out = inflate.flush()
                crc = zlib.crc32(out, crc)
                fh.write(out)
        if part.stat().st_size != member.usize or crc != member.crc:
            part.unlink()
            raise RuntimeError(f"{member.name}: corrupt transfer")
        part.replace(dest)
        return dest


# --- tar streams --------------------------------------------------------------------

def stream_tar(url, handle, progress=None, cancel=None):
    """Read the tar at ``url`` once, front to back, calling ``handle(member, fileobj)``
    for every regular file; ``fileobj`` is the member's content stream."""
    size = content_length(url)
    with _open(url) as resp:
        counted = _Counting(resp, size, progress, cancel)
        with tarfile.open(fileobj=counted, mode="r|") as tar:
            for member in tar:
                if member.isfile():
                    handle(member, tar.extractfile(member))


class _Counting:
    """A read-only file wrapper that reports the bytes consumed."""

    def __init__(self, raw, total, progress, cancel):
        self.raw, self.total, self.progress, self.cancel = raw, total, progress, cancel
        self.done = 0

    def read(self, n=-1):
        if self.cancel is not None:
            self.cancel.check()
        buf = self.raw.read(n)
        self.done += len(buf)
        if self.progress is not None:
            self.progress(self.done, self.total)
        return buf


# --- directory indexes ----------------------------------------------------------------

_ROW_RE = re.compile(r'<a\s+href="([^"?#][^"]*)"[^>]*>.*?</a>(?:.*?<td[^>]*>([^<]*)</td>){2}',
                     re.IGNORECASE | re.DOTALL)
_SIZE_UNITS = {"": 1, "K": 1 << 10, "M": 1 << 20, "G": 1 << 30, "T": 1 << 40}


def _approx_size(text):
    """Bytes from an index page's size column (``182M``, ``6.9K``, ``-`` for a folder)."""
    m = re.fullmatch(r"\s*([\d.]+)\s*([KMGT]?)\s*", text or "")
    return int(float(m.group(1)) * _SIZE_UNITS[m.group(2)]) if m else None


@dataclass(frozen=True)
class IndexEntry:
    name: str
    is_dir: bool
    size: int | None     # rounded as the listing shows it; the transfer reads the exact one


def list_index(url):
    """The entries of an Apache-style directory listing at ``url``, in page order."""
    with _open(url, quick=True) as resp:
        page = resp.read().decode("utf-8", errors="replace")
    out = []
    for row in re.split(r"<tr[ >]", page):
        m = _ROW_RE.search(row)
        if m is None:
            continue
        href = html.unescape(m.group(1))
        if href.startswith("/") or "://" in href:
            continue                # the parent link
        out.append(IndexEntry(href.rstrip("/"), href.endswith("/"), _approx_size(m.group(2))))
    return out
