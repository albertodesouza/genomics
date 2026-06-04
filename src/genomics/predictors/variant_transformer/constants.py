from __future__ import annotations

BASE_TO_ID = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4, "PAD": 5}
ID_TO_BASE = {v: k for k, v in BASE_TO_ID.items()}
BASE_PAD_ID = BASE_TO_ID["PAD"]

VARIANT_TYPE_TO_ID = {"SNP": 0, "INS": 1, "DEL": 2}
ID_TO_VARIANT_TYPE = {v: k for k, v in VARIANT_TYPE_TO_ID.items()}

HAPLOTYPE_TO_ID = {"H1": 0, "H2": 1}
ID_TO_HAPLOTYPE = {v: k for k, v in HAPLOTYPE_TO_ID.items()}

DEFAULT_CLASSES = ["AFR", "AMR", "EAS", "EUR", "SAS"]
