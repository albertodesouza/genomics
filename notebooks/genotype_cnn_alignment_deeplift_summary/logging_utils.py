"""Verbosity control for the shared genotype-based data pipeline."""

from contextlib import contextmanager

from genomics.predictors.genotype_based.data.pipeline import console as _pipeline_console


@contextmanager
def quiet_pipeline_logs():
    """Silences genomics.predictors.genotype_based.data.pipeline's rich-console output (cache
    validation, per-batch recovery progress, etc.) for the duration of the block, without editing
    the shared pipeline module -- that console is a module-level singleton, so this is a
    temporary, restorable flip of its `quiet` flag rather than a monkeypatch of behavior.
    """
    previous = _pipeline_console.quiet
    _pipeline_console.quiet = True
    try:
        yield
    finally:
        _pipeline_console.quiet = previous
