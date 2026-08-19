"""Verbosity control for the shared genotype-based data pipeline."""

from contextlib import contextmanager

from genomics.predictors.genotype_based.data.pipeline import console as _pipeline_console
from genomics.predictors.genotype_based.data.processed_dataset import console as _processed_dataset_console
from genomics.predictors.genotype_based.models.cnn2_model import console as _cnn2_model_console

# Each of these modules instantiates its own module-level rich Console() singleton (not a shared
# one), so silencing pipeline's alone still lets processed_dataset's "Preparando X.pt..." /
# "N samples (lazy)" lines and cnn2_model's "Modelo CNN2 criado: ..." block through.
_QUIETABLE_CONSOLES = (_pipeline_console, _processed_dataset_console, _cnn2_model_console)


@contextmanager
def quiet_pipeline_logs():
    """Silences the rich-console output (cache validation, per-batch recovery progress, model
    construction summaries, etc.) of every module touched by `prepare_data` + CNN2 model
    construction, for the duration of the block, without editing those modules -- each console is
    a module-level singleton, so this is a temporary, restorable flip of each one's `quiet` flag
    rather than a monkeypatch of behavior.
    """
    previous = [c.quiet for c in _QUIETABLE_CONSOLES]
    for c in _QUIETABLE_CONSOLES:
        c.quiet = True
    try:
        yield
    finally:
        for c, was_quiet in zip(_QUIETABLE_CONSOLES, previous):
            c.quiet = was_quiet
