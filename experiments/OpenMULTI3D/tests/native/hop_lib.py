"""The libhop.so that ``make hop_standalone`` copies into the hopacity package, which the
native tests below reach through hopacity. HOP's own tests are in ``HOP/tests``."""
from pathlib import Path

HOP_LIB = Path(__file__).resolve().parents[2] / "HOP" / "hopacity" / "libhop.so"
