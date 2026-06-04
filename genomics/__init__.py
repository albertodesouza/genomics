"""Source-tree compatibility package for the `src/genomics` layout.

This shim lets root-level legacy entrypoints import `genomics.*` directly from a
checkout before the package is installed. Installed environments use the real
package from `src/genomics` via `pyproject.toml`.
"""

from pathlib import Path

_src_package = Path(__file__).resolve().parents[1] / "src" / "genomics"
if _src_package.is_dir():
    __path__.append(str(_src_package))

del Path, _src_package
