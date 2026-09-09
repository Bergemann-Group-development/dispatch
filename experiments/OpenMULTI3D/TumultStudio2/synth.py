"""Compatibility facade for the headless :mod:`tumult.synthesis.workflows`."""
from tumult.synthesis.workflows import *  # noqa: F401,F403
from tumult.synthesis import workflows as _impl


def __getattr__(name):
    return getattr(_impl, name)
