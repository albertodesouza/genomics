#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Legacy wrapper for shared sklearn PCA cache utilities."""

from __future__ import annotations

from genomics.core.sklearn_pca_cache import *  # noqa: F401,F403
from genomics.core.sklearn_pca_cache import run_cli


if __name__ == "__main__":
    run_cli()
