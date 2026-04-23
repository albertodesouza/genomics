from __future__ import annotations

from datasets import Features, Sequence, Value


PREDICTION_MATRIX_FEATURE = Sequence(Sequence(Value("float32")))


def sample_features() -> Features:
    return Features(
        {
            "sample_id": Value("string"),
            "family_id": Value("string"),
            "sex": Value("int32"),
            "population": Value("string"),
            "superpopulation": Value("string"),
            "longevity": Value("bool"),
            "frog_likelihood": Sequence(Value("float32")),
            "frog_population_names": Sequence(Value("string")),
            "available_genes": Sequence(Value("string")),
            "created_at": Value("string"),
            "last_updated": Value("string"),
            "source_dataset": Value("string"),
            "source_path": Value("string"),
            "source_sample_path": Value("string"),
        }
    )


def gene_window_features() -> Features:
    return Features(
        {
            "sample_id": Value("string"),
            "gene": Value("string"),
            "source_dataset": Value("string"),
            "chromosome": Value("string"),
            "start": Value("int64"),
            "end": Value("int64"),
            "window_size": Value("int64"),
            "window_type": Value("string"),
            "ref_sequence": Value("string"),
            "h1_sequence": Value("string"),
            "h2_sequence": Value("string"),
            "outputs": Sequence(Value("string")),
            "ontologies": Sequence(Value("string")),
            "source_path": Value("string"),
            "source_sample_path": Value("string"),
            "source_window_path": Value("string"),
            "source_metadata_path": Value("string"),
            "predictions_h1": {
                "rna_seq": PREDICTION_MATRIX_FEATURE,
                "atac": PREDICTION_MATRIX_FEATURE,
                "cage": PREDICTION_MATRIX_FEATURE,
                "chip_histone": PREDICTION_MATRIX_FEATURE,
                "chip_tf": PREDICTION_MATRIX_FEATURE,
                "dnase": PREDICTION_MATRIX_FEATURE,
                "procap": PREDICTION_MATRIX_FEATURE,
            },
            "predictions_h2": {
                "rna_seq": PREDICTION_MATRIX_FEATURE,
                "atac": PREDICTION_MATRIX_FEATURE,
                "cage": PREDICTION_MATRIX_FEATURE,
                "chip_histone": PREDICTION_MATRIX_FEATURE,
                "chip_tf": PREDICTION_MATRIX_FEATURE,
                "dnase": PREDICTION_MATRIX_FEATURE,
                "procap": PREDICTION_MATRIX_FEATURE,
            },
        }
    )
