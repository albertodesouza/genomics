#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility wrapper for shared sklearn PCA cache utilities.

New code should import from ``genomics_pipeline.sklearn_pca_cache``.
This module remains temporarily so older scripts executed from this directory
continue to work during the migration.
"""

from __future__ import annotations

from genomics_pipeline.sklearn_pca_cache import *  # noqa: F401,F403
from genomics_pipeline.sklearn_pca_cache import run_cli


if __name__ == "__main__":
    run_cli()
