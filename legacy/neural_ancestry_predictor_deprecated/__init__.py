"""Deprecated neural ancestry predictor implementation."""

from __future__ import annotations

import sys


# Historical modules in this package import each other through the old top-level
# package name. Keep that internal alias without restoring root-level wrappers
# or public ``genomics neural`` entrypoints.
sys.modules.setdefault("neural_ancestry_predictor_deprecated", sys.modules[__name__])
