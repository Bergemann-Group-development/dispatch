"""Bokeh directory-app entry point:  bokeh serve --show TumultStudio

Bokeh runs this file as a script with the application directory on ``sys.path``, so
the package itself is reachable only once its parent — the project directory, which
also holds ``tumult`` — is there too. That is the one path bootstrap in the app; no
module inside the package touches ``sys.path``.
"""
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from TumultStudio.app import build  # noqa: E402

build()
