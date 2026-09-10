"""clockcheck — an honest audit harness for epigenetic aging clocks.

Lead with the caveats. Run the checks the papers skip. Never let synthetic
plumbing pass as evidence.
"""
from .clocks import BUILTIN, ClockSpec, LinearClock
from .data import MethylationDataset, load_csv, load_geo, synthetic_fixture
from .report import build_audit_card, card_to_json, card_to_markdown

__version__ = "0.1.0"
__all__ = [
    "LinearClock", "ClockSpec", "BUILTIN",
    "MethylationDataset", "load_csv", "load_geo", "synthetic_fixture",
    "build_audit_card", "card_to_markdown", "card_to_json",
]
