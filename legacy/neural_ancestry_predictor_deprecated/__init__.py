"""Deprecated neural ancestry predictor implementation."""

from __future__ import annotations

import sys


# Historical modules in this package import each other through the old top-level
# package name. Keep that internal alias so ``genomics neural`` works after the
# code moved under ``legacy/``, without restoring root-level wrappers.
sys.modules.setdefault("neural_ancestry_predictor_deprecated", sys.modules[__name__])
