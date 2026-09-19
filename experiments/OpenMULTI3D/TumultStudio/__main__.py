"""Serve TumultStudio behind a landing page:  bin/tpy -m TumultStudio [--port N] [--no-show]

A session's document is built on the server before the first byte of the page is sent,
so a browser opening the app route sees nothing for the whole build. This launcher
serves ``templates/landing.html`` at the site root and opens that instead: the landing
page paints at once and navigates to the app route, which the browser fetches while the
landing page stays on screen. ``bokeh serve --show TumultStudio`` still runs the app
alone, without the landing page.
"""
import argparse
from pathlib import Path

from bokeh.command.util import build_single_handler_application
from bokeh.server.server import Server
from tornado.web import RequestHandler

HERE = Path(__file__).resolve().parent
APP_ROUTE = "/" + HERE.name          # the route bokeh serve gives a directory app


class LandingHandler(RequestHandler):
    """The site root: the static landing page, read on every request so edits show."""

    def get(self):
        self.set_header("Cache-Control", "no-store")
        self.write((HERE / "templates" / "landing.html").read_text())


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m TumultStudio",
                                     description="Serve TumultStudio behind a landing page")
    parser.add_argument("--port", type=int, default=5006)
    parser.add_argument("--no-show", action="store_true",
                        help="do not open the landing page in a browser")
    args = parser.parse_args(argv)
    application = build_single_handler_application(str(HERE))
    # no bokeh index at "/": the landing page takes the root
    server = Server({APP_ROUTE: application}, port=args.port, use_index=False,
                    extra_patterns=[("/", LandingHandler)])
    server.start()
    if not args.no_show:
        server.io_loop.add_callback(server.show, "/")
    print(f"TumultStudio at http://localhost:{args.port}/")
    server.io_loop.start()


if __name__ == "__main__":
    main()
