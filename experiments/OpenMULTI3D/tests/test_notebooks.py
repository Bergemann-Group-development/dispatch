"""Guard the de-facto notebook API.

The demo notebooks are the library's public contract, but executing them shells out to
``./dispatch.x`` (minutes of Fortran) - so the default test statically extracts every
``tumult.<name>`` / ``tm.<name>`` / ``from tumult... import`` usage from the notebooks
and asserts each name still resolves. This is the tripwire for Phase 1's rewrite of
tumult/__init__.py (wildcard imports -> explicit __all__).

Full execution is opt-in:  TUMULT_EXECUTE_NOTEBOOKS=1 pytest tests/test_notebooks.py
"""
import ast
import json
import os
import subprocess
import sys

import pytest

from conftest import ROOT

NOTEBOOKS = [ROOT / "demo.ipynb", ROOT / "demo_corrections.ipynb", ROOT / "hop.ipynb"]


def _code_cells(path):
    cells = json.loads(path.read_text())["cells"]
    for cell in cells:
        if cell["cell_type"] == "code":
            src = "".join(cell["source"])
            # drop IPython magics / shell escapes that ast can't parse
            yield "\n".join(l for l in src.splitlines()
                            if not l.lstrip().startswith(("!", "%")))


def _tumult_names_used(paths):
    """All attribute names accessed on the tumult package (or its `tm` alias) plus all
    `from tumult[.mod] import name` imports across the notebooks."""
    attrs, imports = set(), set()
    for path in paths:
        for src in _code_cells(path):
            try:
                tree = ast.parse(src)
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
                        and node.value.id in ("tumult", "tm")):
                    attrs.add(node.attr)
                elif isinstance(node, ast.ImportFrom) and node.module \
                        and node.module.split(".")[0] == "tumult":
                    for alias in node.names:
                        imports.add((node.module, alias.name))
    return attrs, imports


@pytest.mark.parametrize("nb", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_exists(nb):
    assert nb.exists()


def test_notebook_api_surface_resolves():
    import importlib

    import tumult

    attrs, imports = _tumult_names_used(NOTEBOOKS)
    assert attrs or imports, "extraction found nothing - parser broken?"

    missing = sorted(a for a in attrs if not hasattr(tumult, a))
    assert not missing, (
        f"notebook-used names no longer resolve on the tumult package: {missing}"
    )
    for module, name in sorted(imports):
        mod = importlib.import_module(module)
        assert name == "*" or hasattr(mod, name), f"from {module} import {name} broke"


@pytest.mark.skipif(not os.environ.get("TUMULT_EXECUTE_NOTEBOOKS"),
                    reason="set TUMULT_EXECUTE_NOTEBOOKS=1 to execute (runs dispatch.x)")
@pytest.mark.parametrize("nb", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_executes(nb):
    subprocess.run(
        [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
         "--stdout", str(nb)],
        cwd=ROOT, check=True, stdout=subprocess.DEVNULL, timeout=3600,
    )
