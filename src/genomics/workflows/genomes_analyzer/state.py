"""Shared runtime state for the legacy pipeline implementation.

The heavy pipeline still stores mutable process-wide state in the legacy module.
These helpers keep the public package API stable while the internals are split
incrementally.
"""

from __future__ import annotations

from . import legacy


def set_config(config: dict) -> None:
    legacy.cfg_global = config


def get_config() -> dict | None:
    return legacy.cfg_global


def get_console():
    return legacy.console
