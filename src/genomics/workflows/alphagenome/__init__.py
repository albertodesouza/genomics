"""AlphaGenome analysis workflow.

The analyzer names are imported lazily (PEP 562): ``neural_module`` needs rich and the AlphaGenome
client, which light submodules such as ``outputs`` and ``catalog`` (used by the visualizer) do not.
"""
from importlib import import_module
from typing import Any

__all__ = ["AlphaGenomeAnalyzer", "DEFAULT_CONFIG", "parse_fasta", "validate_sequence"]


def __getattr__(name: str) -> Any:
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(".neural_module", __name__), name)
    globals()[name] = value
    return value
