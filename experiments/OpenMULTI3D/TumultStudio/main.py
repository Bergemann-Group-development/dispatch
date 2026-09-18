"""Bokeh directory-app entry point:  bokeh serve --show TumultStudio

Bokeh runs this file as a script with the application directory on ``sys.path``, so
the package itself is reachable only once its parent — the project directory, which
also holds ``tumult`` — is there too. This and ``app_hooks.py`` are the only path
bootstraps in the app; no module inside the package touches ``sys.path``. Bokeh restores
``sys.path`` once this script has run, which is harmless here: ``build()`` runs inside
the script, and every later import in the package is of a ``tumult`` or ``TumultStudio``
submodule, resolved through the already-imported package's ``__path__``.
"""
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from TumultStudio.app import build  # noqa: E402

build()
