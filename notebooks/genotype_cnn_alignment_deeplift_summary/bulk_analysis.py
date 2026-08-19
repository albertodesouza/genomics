"""Population-scale (offline, pre-computed) knockout results loaded from CSV."""

import pandas as pd


def load_bulk_results(bulk_results_path):
    """Loads the offline bulk knockout CSV, derives log-odds columns, and excludes the disabled
    alphagenome_atac method (see the notebook's Ontology caveat markdown) from every downstream
    stat/plot."""
    bulk_df = pd.read_csv(bulk_results_path)
    bulk_df = bulk_df[bulk_df["method"] != "alphagenome_atac"]
    # log-odds(strong/weak) = logit_strong - logit_weak = log(p_strong / p_weak) -- exact (not an
    # approximation) for this binary softmax head, since the normalizing constant cancels in the
    # difference. The bulk CSV only stores the two raw per-class logits, so it's derived here
    # rather than re-running the (separate, offline) knockout script that produced the CSV.
    bulk_df["baseline_log_odds"] = bulk_df["baseline_strong_logit"] - bulk_df["baseline_weak_logit"]
    bulk_df["perturbed_log_odds"] = bulk_df["perturbed_strong_logit"] - bulk_df["perturbed_weak_logit"]
    bulk_df["delta_log_odds"] = bulk_df["perturbed_log_odds"] - bulk_df["baseline_log_odds"]
    return bulk_df


def compute_gene_method_stats(bulk_df):
    gene_method_stats = (
        bulk_df.groupby(["gene", "method"])
        .agg(
            mean_baseline=("baseline_log_odds", "mean"),
            mean_perturbed=("perturbed_log_odds", "mean"),
            mean_delta=("delta_log_odds", "mean"),
            std_delta=("delta_log_odds", "std"),
            flip_rate=("flipped", "mean"),
            n=("flipped", "size"),
        )
        .reset_index()
    )
    gene_method_stats["flip_rate_pct"] = gene_method_stats["flip_rate"] * 100
    return gene_method_stats
