"""All matplotlib figure builders used by the notebook."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.metrics import confusion_matrix

from .annotations import (
    gene_boundary_positions, gene_exon_haplotype_local_boxes, gene_exon_local_boxes,
    get_gene_tss, window_genomic_axis,
)
from .attribution import class_track_curves_mean_std, gene_track_curves, smooth
from .cage import find_cage_promoter_individual, find_cage_promoter_no_tss_individual
from .knockout import cnn_crop_bounds


def plot_confusion_matrices(predictions, class_names, title, ncols, show_accuracy=False):
    """Confusion-matrix grid, one panel per model in `predictions` ({name: {"y_true", "y_pred"}}).
    `ncols` controls the grid width (nrows is the ceiling of n_models / ncols); pass
    show_accuracy=True to include each panel's accuracy in its title."""
    names = list(predictions)
    n_models = len(names)
    nrows = -(-n_models // ncols)  # ceil
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.5 * nrows))
    axes = np.atleast_1d(axes).flatten()
    for ax, name in zip(axes, names):
        y_true, y_pred = predictions[name]["y_true"], predictions[name]["y_pred"]
        cm = confusion_matrix(y_true, y_pred)
        ax.imshow(cm, cmap="Blues")
        label = name
        if show_accuracy:
            acc = (y_true == y_pred).mean()
            label = f"{name} (acc={acc:.3f})"
        ax.set_title(label)
        ax.set_xticks(range(len(class_names)))
        ax.set_xticklabels(class_names, rotation=30, ha="right")
        ax.set_yticks(range(len(class_names)))
        ax.set_yticklabels(class_names)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                        color="white" if cm[i, j] > cm.max() / 2 else "black")
    for ax in axes[n_models:]:
        ax.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


def plot_signal_std_bands(genes, config, target_idx, other_idx, class_names, class_colors,
                           ontology_terms, strands, dataset_dir, window_center_size,
                           mean_input_total, loaders, model_name, gene_id_map,
                           transcript_extractor_mane, transcript_extractor_coding, out_dir):
    """Per-gene RNA-seq signal (mean +/-1 SD across individuals) for the target vs. other class,
    one figure per gene, saved to `out_dir` and returned as (figs, labels) for a carousel."""
    out_dir.mkdir(parents=True, exist_ok=True)
    n_ont = len(ontology_terms)
    n_rows = n_ont + 1  # RNA-seq signal row per ontology + exon row (no DeepLIFT row here)

    figs, labels = [], []
    for gene in genes:
        chrom, window_start_0based, local_length = window_genomic_axis(dataset_dir, gene, window_center_size)
        exon_boxes, gene_strand = gene_exon_local_boxes(
            gene, chrom, window_start_0based, local_length,
            gene_id_map, transcript_extractor_mane, transcript_extractor_coding,
        )
        tss_local, tes_local = gene_boundary_positions(exon_boxes, gene_strand)

        mean_std_by_class = {
            target_idx: class_track_curves_mean_std(loaders, model_name, target_idx, gene, config),
            other_idx: class_track_curves_mean_std(loaders, model_name, other_idx, gene, config),
        }
        total_signal_curves = gene_track_curves(mean_input_total, config, gene)

        fig, axes = plt.subplots(
            n_rows, 2, figsize=(13, 1.6 * n_rows + 1.0), sharex=True,
            gridspec_kw={"height_ratios": [2] * n_ont + [1]},
        )
        for row in range(n_rows - 1):
            axes[row, 0].sharey(axes[row, 1])

        for col, strand in enumerate(strands):
            is_sense = strand == gene_strand
            for i, ontology in enumerate(ontology_terms):
                ax = axes[i, col]
                for class_idx in (target_idx, other_idx):
                    pair = mean_std_by_class[class_idx].get((ontology, strand))
                    if pair is None:
                        continue
                    mean_curve, std_curve = pair
                    mean_s, std_s = smooth(mean_curve), smooth(std_curve)
                    x = np.arange(len(mean_s))
                    color = class_colors[class_idx]
                    ax.plot(x, mean_s, color=color, linewidth=1.1, label=class_names[class_idx])
                    # Solid, semi-transparent fill instead of a hatched outline -- where the two
                    # classes' bands overlap, matplotlib's alpha compositing blends red and blue
                    # into a visibly distinct mixed color, so overlap reads directly from the
                    # background instead of needing to spot two crossing hatch patterns.
                    ax.fill_between(x, mean_s - std_s, mean_s + std_s, color=color, alpha=0.3,
                                     linewidth=0.0, zorder=0)
                total_curve = total_signal_curves.get((ontology, strand))
                if total_curve is not None:
                    ax.plot(smooth(total_curve), color="#808080", linewidth=0.8, linestyle="--", label="all individuals")
                if is_sense:
                    for s, e in exon_boxes:
                        ax.axvspan(s, e, color="tab:green", alpha=0.10, zorder=0)
                    if tss_local is not None:
                        ax.axvline(tss_local, color="tab:green", linewidth=1.0, linestyle=":", alpha=0.8)
                        ax.axvline(tes_local, color="tab:green", linewidth=1.0, linestyle="-.", alpha=0.8)
                        ax.text(tss_local, 1.02, "TSS", transform=ax.get_xaxis_transform(),
                                fontsize=6, color="tab:green", ha="center", va="bottom")
                        ax.text(tes_local, 1.02, "TES", transform=ax.get_xaxis_transform(),
                                fontsize=6, color="tab:green", ha="center", va="bottom")
                if col == 0:
                    ax.set_ylabel(f"{ontology}\nRNA-seq signal\n(normalized)", fontsize=7)
                if i == 0:
                    suffix = " (sense)" if is_sense else " (antisense)"
                    ax.set_title(f"{strand} strand{suffix}", fontsize=9)
                    if col == 0:
                        ax.legend(fontsize=6, loc="upper right")

            ax_exon = axes[n_rows - 1, col]
            if is_sense:
                for s, e in exon_boxes:
                    ax_exon.add_patch(plt.Rectangle((s, 0.15), e - s, 0.7, color="tab:green", alpha=0.8))
                ax_exon.set_ylabel(f"exons\n({strand})", fontsize=7, rotation=0, labelpad=22, va="center")
            else:
                ax_exon.text(0.5, 0.5, "exons not applicable\n(antisense)", transform=ax_exon.transAxes,
                              ha="center", va="center", fontsize=7, color="gray")
            ax_exon.set_xlim(0, local_length)
            ax_exon.set_ylim(0, 1)
            ax_exon.set_yticks([])
            ax_exon.set_xlabel(f"Position in window (bp) -- start={chrom}:{window_start_0based + 1:,}", fontsize=7.5)

        fig.suptitle(
            f"{gene} -- RNA-seq signal: {class_names[target_idx]} vs. {class_names[other_idx]} "
            f"(shaded = mean +/-1 SD across individuals, mixed color = overlap), gene strand={gene_strand}"
        )
        fig.tight_layout()
        out_path = out_dir / f"signal_std_band_{gene}.png"
        fig.savefig(out_path, dpi=140)
        plt.close(fig)
        figs.append(fig)
        labels.append(gene)
        print(f"saved {out_path}")

    return figs, labels


def plot_exon_overlap(genes, config, target_idx, other_idx, class_names, class_colors,
                       ontology_terms, strands, dataset_dir, window_center_size,
                       deeplift_results, mean_input_total, logodds_attr_target_pop,
                       logodds_attr_other_pop, gene_id_map, transcript_extractor_mane,
                       transcript_extractor_coding, out_dir):
    """Per-gene (RNA-seq signal, DeepLIFT log-odds attribution) row pairs plus an exon track, one
    figure per gene, saved to `out_dir` and returned as (figs, labels) for a carousel. Also
    returns the per-track Pearson r between the smoothed strong-minus-weak signal diff and the
    smoothed log-odds attribution curve, as a DataFrame."""
    out_dir.mkdir(parents=True, exist_ok=True)
    n_ont = len(ontology_terms)
    n_rows = 2 * n_ont + 1  # (RNA-seq signal, DeepLIFT attribution) pair per ontology, + exon row

    figs, labels = [], []
    all_corr_rows = []
    for gene in genes:
        chrom, window_start_0based, local_length = window_genomic_axis(dataset_dir, gene, window_center_size)
        exon_boxes, gene_strand = gene_exon_local_boxes(
            gene, chrom, window_start_0based, local_length,
            gene_id_map, transcript_extractor_mane, transcript_extractor_coding,
        )
        tss_local, tes_local = gene_boundary_positions(exon_boxes, gene_strand)

        # One curve per population: how much that population's genotype pushes the model's
        # log-odds(target/other) up (toward target) or down (toward other) at each window position.
        target_pop_logodds_curves = gene_track_curves(logodds_attr_target_pop, config, gene)
        other_pop_logodds_curves = gene_track_curves(logodds_attr_other_pop, config, gene)

        class_signal_curves = {
            target_idx: gene_track_curves(deeplift_results[target_idx]["mean_input"], config, gene),
            other_idx: gene_track_curves(deeplift_results[other_idx]["mean_input"], config, gene),
        }
        total_signal_curves = gene_track_curves(mean_input_total, config, gene)

        fig, axes = plt.subplots(
            n_rows, 2, figsize=(13, 1.5 * n_rows + 1.0), sharex=True,
            gridspec_kw={"height_ratios": [2, 2] * n_ont + [1]},
        )
        for row in range(n_rows - 1):
            axes[row, 0].sharey(axes[row, 1])

        corr_rows = []
        for col, strand in enumerate(strands):
            is_sense = strand == gene_strand
            for i, ontology in enumerate(ontology_terms):
                signal_row, attr_row = 2 * i, 2 * i + 1

                ax_sig = axes[signal_row, col]
                target_curve = class_signal_curves[target_idx].get((ontology, strand))
                other_curve = class_signal_curves[other_idx].get((ontology, strand))
                total_curve = total_signal_curves.get((ontology, strand))
                diff_smoothed = None
                ax_diff = None
                if target_curve is not None and other_curve is not None:
                    target_s, other_s = smooth(target_curve), smooth(other_curve)
                    x = np.arange(len(target_s))
                    # Shade the gap between the two class means so it reads at a glance even
                    # where the raw curves nearly overlap -- dark where target sits above other,
                    # light where it's the other way around.
                    ax_sig.fill_between(x, target_s, other_s, where=target_s >= other_s,
                                         color=class_colors[target_idx], alpha=0.15, interpolate=True, zorder=0)
                    ax_sig.fill_between(x, target_s, other_s, where=target_s < other_s,
                                         color=class_colors[other_idx], alpha=0.15, interpolate=True, zorder=0)
                    ax_sig.plot(x, target_s, color=class_colors[target_idx], linewidth=1.0, label=class_names[target_idx])
                    ax_sig.plot(x, other_s, color=class_colors[other_idx], linewidth=1.0, label=class_names[other_idx])
                    diff_smoothed = target_s - other_s

                    # The diff is typically tiny next to the raw signal, so sharing the main axis
                    # squashes it flat and makes its sign unreadable. Give it its own (per-column,
                    # unshared) y-scale via a twin axis instead, anchored by a dashed zero-line.
                    ax_diff = ax_sig.twinx()
                    ax_diff.axhline(0, color="tab:purple", linewidth=0.6, linestyle="--", alpha=0.5, zorder=1)
                    ax_diff.plot(x, diff_smoothed, color="tab:purple", linewidth=1.0, linestyle=":",
                                 label=f"{class_names[target_idx]} - {class_names[other_idx]} (diff)")
                    ax_diff.tick_params(axis="y", labelsize=6, labelcolor="tab:purple")
                    if col == 1:
                        ax_diff.set_ylabel("diff\n(target - other)", fontsize=7, color="tab:purple")
                if total_curve is not None:
                    ax_sig.plot(smooth(total_curve), color="#808080", linewidth=0.8, linestyle="--",
                                label="all individuals", zorder=1)
                if is_sense:
                    for s, e in exon_boxes:
                        ax_sig.axvspan(s, e, color="tab:green", alpha=0.10, zorder=0)
                    if tss_local is not None:
                        ax_sig.axvline(tss_local, color="tab:green", linewidth=1.0, linestyle=":", alpha=0.8)
                        ax_sig.axvline(tes_local, color="tab:green", linewidth=1.0, linestyle="-.", alpha=0.8)
                        ax_sig.text(tss_local, 1.02, "TSS", transform=ax_sig.get_xaxis_transform(),
                                    fontsize=6, color="tab:green", ha="center", va="bottom")
                        ax_sig.text(tes_local, 1.02, "TES", transform=ax_sig.get_xaxis_transform(),
                                    fontsize=6, color="tab:green", ha="center", va="bottom")
                if col == 0:
                    ax_sig.set_ylabel(f"{ontology}\nRNA-seq signal\n(normalized)", fontsize=7)
                if signal_row == 0:
                    suffix = " (sense)" if is_sense else " (antisense)"
                    ax_sig.set_title(f"{strand} strand{suffix}", fontsize=9)
                    if col == 0:
                        handles, plot_labels = ax_sig.get_legend_handles_labels()
                        if ax_diff is not None:
                            diff_handles, diff_labels = ax_diff.get_legend_handles_labels()
                            handles, plot_labels = handles + diff_handles, plot_labels + diff_labels
                        ax_sig.legend(handles, plot_labels, fontsize=6, loc="upper right")

                # Both curves share one axis and one sign convention: positive = pushes the
                # model's log-odds(target/other) toward target, negative = toward other,
                # regardless of which population is plotted. Color identifies the population.
                ax_attr = axes[attr_row, col]
                target_pop_curve = target_pop_logodds_curves.get((ontology, strand))
                if target_pop_curve is not None:
                    ax_attr.plot(smooth(target_pop_curve), color=class_colors[target_idx],
                                 linewidth=1.0, linestyle="-",
                                 label=f"{class_names[target_idx]} individuals")
                other_pop_curve = other_pop_logodds_curves.get((ontology, strand))
                if other_pop_curve is not None:
                    ax_attr.plot(smooth(other_pop_curve), color=class_colors[other_idx],
                                 linewidth=1.0, linestyle="-",
                                 label=f"{class_names[other_idx]} individuals")
                ax_attr.axhline(0, color="black", linewidth=0.7, alpha=0.5, zorder=1)
                if is_sense:
                    for s, e in exon_boxes:
                        ax_attr.axvspan(s, e, color="tab:green", alpha=0.10, zorder=0)
                    if tss_local is not None:
                        ax_attr.axvline(tss_local, color="tab:green", linewidth=1.0, linestyle=":", alpha=0.8)
                        ax_attr.axvline(tes_local, color="tab:green", linewidth=1.0, linestyle="-.", alpha=0.8)
                if col == 0:
                    ax_attr.set_ylabel(f"{ontology}\nDeepLIFT attr\ntoward log-odds", fontsize=7)
                    if i == 0:
                        ax_attr.legend(fontsize=6, loc="upper right")

                # Quantify -- rather than just eyeball -- whether where the two classes' RNA-seq
                # signal differs lines up with where each population's log-odds attribution is
                # large: Pearson r between the smoothed diff and the smoothed attribution curve.
                if diff_smoothed is not None:
                    if other_pop_curve is not None:
                        r = float(np.corrcoef(diff_smoothed, smooth(other_pop_curve))[0, 1])
                        corr_rows.append({
                            "strand": strand, "ontology": ontology, "individuals": class_names[other_idx],
                            "corr(diff_signal, logodds_attribution)": r,
                        })
                    if target_pop_curve is not None:
                        r = float(np.corrcoef(diff_smoothed, smooth(target_pop_curve))[0, 1])
                        corr_rows.append({
                            "strand": strand, "ontology": ontology, "individuals": class_names[target_idx],
                            "corr(diff_signal, logodds_attribution)": r,
                        })

            ax_exon = axes[n_rows - 1, col]
            if is_sense:
                for s, e in exon_boxes:
                    ax_exon.add_patch(plt.Rectangle((s, 0.15), e - s, 0.7, color="tab:green", alpha=0.8))
                ax_exon.set_ylabel(f"exons\n({strand})", fontsize=7, rotation=0, labelpad=22, va="center")
            else:
                ax_exon.text(0.5, 0.5, "exons not applicable\n(antisense)", transform=ax_exon.transAxes,
                              ha="center", va="center", fontsize=7, color="gray")
            ax_exon.set_xlim(0, local_length)
            ax_exon.set_ylim(0, 1)
            ax_exon.set_yticks([])
            ax_exon.set_xlabel(f"Position in window (bp) -- start={chrom}:{window_start_0based + 1:,}", fontsize=7.5)

        fig.suptitle(
            f"{gene} -- gene strand={gene_strand} ({len(exon_boxes)} exons), "
            f"positive = pushes the log-odds toward {class_names[target_idx]}"
        )
        fig.tight_layout()
        out_path = out_dir / f"exon_overlap_{gene}.png"
        fig.savefig(out_path, dpi=140)
        plt.close(fig)
        figs.append(fig)
        labels.append(gene)
        all_corr_rows.extend({"gene": gene, **row} for row in corr_rows)
        print(f"saved {out_path}")

    return figs, labels, pd.DataFrame(all_corr_rows)


def plot_cage_gene_start_centered_individual(ctx, sample_id, gene):
    _, _, strand = get_gene_tss(gene, ctx.tss_df_mane, ctx.tss_df_coding)
    dark_color = ctx.class_colors[ctx.target_idx]
    light_color = ctx.class_colors[ctx.other_idx]
    fig, axes = plt.subplots(
        2, 2, figsize=(13, 5.5), sharey="row", gridspec_kw={"height_ratios": [3, 1]},
    )
    for col, hap in enumerate(("H1", "H2")):
        tss_result = find_cage_promoter_individual(ctx, sample_id, gene, hap)
        gene_start_result = find_cage_promoter_no_tss_individual(ctx, sample_id, gene, hap)
        gene_start_anchor = gene_start_result["anchor_local_idx"]
        tss_anchor = tss_result["anchor_local_idx"]

        # Display window: the union of both methods' own +/-5 kb search windows, so both
        # neighborhoods stay visible regardless of how far apart the TSS and gene-start anchors
        # sit (minus-strand genes can put these tens of kb apart).
        lo = min(tss_result["search_lo"], gene_start_result["search_lo"])
        hi = max(tss_result["search_hi"], gene_start_result["search_hi"])

        exon_boxes, exon_strand = gene_exon_haplotype_local_boxes(
            ctx.dataset_dir, sample_id, gene, hap, ctx.gene_id_map,
            ctx.transcript_extractor_mane, ctx.transcript_extractor_coding,
        )
        exon_tss_local, exon_tes_local = gene_boundary_positions(exon_boxes, exon_strand)

        ax = axes[0, col]
        x = np.arange(lo, hi) - gene_start_anchor
        # log1p, not raw signal -- same convention the RNA-seq input tracks already use, and
        # standard for this kind of long-tailed, mostly-near-zero signal: keeps zero-valued
        # positions (most of the window, outside real initiation sites) well-defined while
        # compressing the CAGE summit's large dynamic range so both are visible on one axis.
        ax.plot(x, np.log1p(gene_start_result["dark_curve"][lo:hi]), color=dark_color,
                linewidth=1.1, label="dark melanocyte")
        ax.plot(x, np.log1p(gene_start_result["light_curve"][lo:hi]), color=light_color,
                linewidth=1.1, label="light melanocyte")

        ax.axvline(0, color="tab:purple", linestyle="-", linewidth=1.2, alpha=0.6, label="gene start (center)")
        ax.axvline(tss_anchor - gene_start_anchor, color="tab:blue", linestyle="-", alpha=0.6, linewidth=1.4,
                   label=f"TSS ({strand} strand)")
        ax.axvline(tss_result["summit_local_idx"] - gene_start_anchor, color="tab:green",
                   linestyle="--", linewidth=1.4, alpha=0.6, label="CAGE summit (TSS neighborhood)")
        ax.axvline(gene_start_result["summit_local_idx"] - gene_start_anchor, color="tab:orange",
                   linestyle=":", linewidth=1.8, alpha=0.6, label="CAGE summit (gene-start neighborhood)")
        if exon_boxes:
            for s, e in exon_boxes:
                ax.axvspan(s - gene_start_anchor, e - gene_start_anchor, color="tab:green", alpha=0.10, zorder=0)
            if exon_tss_local is not None:
                ax.axvline(exon_tss_local - gene_start_anchor, color="tab:green", linewidth=1.0,
                           linestyle=":", alpha=0.8)
                ax.axvline(exon_tes_local - gene_start_anchor, color="tab:green", linewidth=1.0,
                           linestyle="-.", alpha=0.8)
                ax.text(exon_tss_local - gene_start_anchor, 1.02, "TSS", transform=ax.get_xaxis_transform(),
                        fontsize=6, color="tab:green", ha="center", va="bottom")
                ax.text(exon_tes_local - gene_start_anchor, 1.02, "TES", transform=ax.get_xaxis_transform(),
                        fontsize=6, color="tab:green", ha="center", va="bottom")

        ax.set_title(hap)
        if col == 0:
            ax.set_ylabel("ln(1 + CAGE signal)")
            ax.legend(fontsize=7, loc="upper right")

        ax_exon = axes[1, col]
        if exon_boxes:
            for s, e in exon_boxes:
                ax_exon.add_patch(plt.Rectangle((s - gene_start_anchor, 0.15), e - s, 0.7,
                                                 color="tab:green", alpha=0.8))
            ax_exon.set_ylabel(f"exons\n({exon_strand})", fontsize=7, rotation=0, labelpad=22, va="center")
        else:
            ax_exon.text(0.5, 0.5, "no exons found", transform=ax_exon.transAxes,
                          ha="center", va="center", fontsize=7, color="gray")
        ax_exon.set_xlim(x[0], x[-1])
        ax_exon.set_ylim(0, 1)
        ax_exon.set_yticks([])
        ax_exon.set_xlabel("Position relative to gene start (bp)")

    fig.suptitle(f"{gene} -- {sample_id} -- individual CAGE signal centered on gene start")
    fig.tight_layout()
    out_path = ctx.fig_dir / f"cage_gene_start_centered_{gene}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return fig


def cage_tss_distance_table(ctx, sample_id, genes):
    """Distance (bp) between each technique's candidate promoter location and the literature TSS
    (GENCODE MANE Select), per gene/haplotype: the strand-naive gene start, the CAGE summit found
    in the TSS neighborhood, and the CAGE summit found in the gene-start neighborhood."""
    rows = []
    for gene in genes:
        for hap in ("H1", "H2"):
            tss_result = find_cage_promoter_individual(ctx, sample_id, gene, hap)
            gene_start_result = find_cage_promoter_no_tss_individual(ctx, sample_id, gene, hap)
            tss_anchor = tss_result["anchor_local_idx"]
            rows.extend([
                {"gene": gene, "haplotype": hap, "technique": "gene_start (anchor)",
                 "distance_from_tss_bp": gene_start_result["anchor_local_idx"] - tss_anchor},
                {"gene": gene, "haplotype": hap, "technique": "cage_summit_tss_neighborhood",
                 "distance_from_tss_bp": tss_result["summit_local_idx"] - tss_anchor},
                {"gene": gene, "haplotype": hap, "technique": "cage_summit_gene_start_neighborhood",
                 "distance_from_tss_bp": gene_start_result["summit_local_idx"] - tss_anchor},
            ])
    return pd.DataFrame(rows)


def plot_knockout_result(ctx, result, marker_margin=1000):
    # Each of the 6 tracks (3 ontology terms x 2 strands) plotted separately, not averaged --
    # they come from different cell types and strands, so a mean across them isn't a meaningful
    # biological quantity (a strong effect on one track can be washed out by four flat ones).
    track_meta = result["hap_tracks"]["H1"]["meta"]
    track_labels = [f"{m['ontology_curie']} ({m['strand']})" for m in track_meta]
    n_tracks = len(track_labels)

    orig_len = result["hap_tracks"]["H1"]["orig_array"].shape[0]
    crop_lo, crop_hi = cnn_crop_bounds(orig_len, ctx.config.dataset_input.window_center_size)

    # sharey="row": H1 vs H2 for the *same* track must be on the same y-scale to be visually
    # comparable -- without it matplotlib auto-scales each subplot independently, so a flat H2
    # track can look deceptively different from an active H1 track just because of axis scaling.
    # One extra row (n_tracks + 1) for the exon track, same convention as plot_exon_overlap.
    fig, axes = plt.subplots(
        n_tracks + 1, 2, figsize=(13, 1.9 * n_tracks + 1.3), sharex=True, sharey="row",
        gridspec_kw={"height_ratios": [1] * n_tracks + [0.4]},
    )
    for col, hap in enumerate(["H1", "H2"]):
        track = result["hap_tracks"][hap]
        orig, mod = track["orig_array"], track["mod_array"]
        # Anchor x=0 at *this haplotype's own* indel-corrected TSS -- H1 and H2 can carry
        # different indels upstream of the target, so "the TSS" can sit tens to hundreds of bp
        # apart between the two haplotypes' actual consensus sequences.
        hap_tss_idx = track["marker_local_idx"]["biology_tss"]
        # Base window = the CNN's actual 32,768 bp input crop (so what's plotted matches what the
        # model sees) -- extended only as needed to include this haplotype's own TSS and/or CAGE
        # summit (each with a 1 kb margin) when a marker falls outside that crop. Computed per
        # haplotype now (not from a shared reference-coordinate estimate) since cage_melanocyte's
        # summit is found natively per haplotype and can genuinely differ between H1 and H2.
        hap_markers = [hap_tss_idx, track["marker_local_idx"]["cage_melanocyte"]]
        win_lo = min([crop_lo] + [m - marker_margin for m in hap_markers if m < crop_lo])
        win_hi = max([crop_hi] + [m + marker_margin for m in hap_markers if m >= crop_hi])
        lo, hi = max(0, win_lo), min(orig.shape[0], win_hi)
        x = np.arange(lo, hi) - hap_tss_idx

        exon_boxes, exon_strand = gene_exon_haplotype_local_boxes(
            ctx.dataset_dir, result["sample_id"], result["gene"], hap, ctx.gene_id_map,
            ctx.transcript_extractor_mane, ctx.transcript_extractor_coding,
        )
        exon_tss_local, exon_tes_local = gene_boundary_positions(exon_boxes, exon_strand)

        for row in range(n_tracks):
            ax = axes[row, col]
            # Shaded band = the CNN's actual 32 kb input crop, so it's visually obvious how much
            # of the plotted window (if any extension was needed) the model never directly sees.
            ax.axvspan(crop_lo - hap_tss_idx, crop_hi - hap_tss_idx, color="gray", alpha=0.10, zorder=0)
            ax.plot(x, orig[lo:hi, row], color="tab:blue", linewidth=0.6, alpha=0.65, label="original")
            ax.plot(x, mod[lo:hi, row], color="tab:red", linewidth=0.6, alpha=0.65, label="knockout")
            # Both candidate locations, on every panel, at *this haplotype's* own position --
            # solid+thicker for whichever one this particular experiment actually edited (plus a
            # shaded band showing the *actual* scrambled window for that method), thin+dotted for
            # the other method's site.
            marker_colors = {"biology_tss": "#1f77b4", "cage_melanocyte": "#2ca02c"}
            for marker_method, marker_color in marker_colors.items():
                marker_idx = track["marker_local_idx"][marker_method]
                active = marker_method == result["method"]
                if active:
                    ax.axvspan(track["scramble_start"] - hap_tss_idx, track["scramble_end"] - hap_tss_idx,
                               color=marker_color, alpha=0.20, zorder=1)
                    ax.axvline(marker_idx - hap_tss_idx, color=marker_color, linewidth=1.4, alpha=0.9)
                else:
                    ax.axvline(marker_idx - hap_tss_idx, color=marker_color, linewidth=0.9, linestyle=":", alpha=0.55)
            if exon_boxes:
                for s, e in exon_boxes:
                    ax.axvspan(s - hap_tss_idx, e - hap_tss_idx, color="tab:green", alpha=0.10, zorder=0)
                if exon_tss_local is not None:
                    ax.axvline(exon_tss_local - hap_tss_idx, color="tab:green", linewidth=1.0,
                               linestyle=":", alpha=0.8)
                    ax.axvline(exon_tes_local - hap_tss_idx, color="tab:green", linewidth=1.0,
                               linestyle="-.", alpha=0.8)
                    if row == 0:
                        ax.text(exon_tss_local - hap_tss_idx, 1.02, "TSS", transform=ax.get_xaxis_transform(),
                                fontsize=6, color="tab:green", ha="center", va="bottom")
                        ax.text(exon_tes_local - hap_tss_idx, 1.02, "TES", transform=ax.get_xaxis_transform(),
                                fontsize=6, color="tab:green", ha="center", va="bottom")
            if col == 0:
                ax.set_ylabel(track_labels[row], fontsize=7.5)
            if row == 0:
                ax.set_title(
                    f"{hap} (scrambled {track['scramble_start'] - hap_tss_idx}..{track['scramble_end'] - hap_tss_idx})"
                )
                handles = [
                    plt.Line2D([], [], color="tab:blue", linewidth=0.6, alpha=0.65, label="original"),
                    plt.Line2D([], [], color="tab:red", linewidth=0.6, alpha=0.65, label="knockout"),
                    plt.Line2D([], [], color="#1f77b4", linewidth=1.6, label="TSS"),
                    plt.Line2D([], [], color="#2ca02c", linewidth=1.6, label="CAGE summit"),
                    Patch(facecolor="gray", alpha=0.10, label="CNN 32 kb input crop"),
                    Patch(facecolor="gray", alpha=0.20, label="scrambled window (active method)"),
                    Patch(facecolor="tab:green", alpha=0.10, label="exon"),
                ]
                ax.legend(handles=handles, fontsize=5.5, ncol=2)

        ax_exon = axes[n_tracks, col]
        if exon_boxes:
            for s, e in exon_boxes:
                ax_exon.add_patch(plt.Rectangle((s - hap_tss_idx, 0.15), e - s, 0.7, color="tab:green", alpha=0.8))
            ax_exon.set_ylabel(f"exons\n({exon_strand})", fontsize=7, rotation=0, labelpad=22, va="center")
        else:
            ax_exon.text(0.5, 0.5, "no exons found", transform=ax_exon.transAxes,
                          ha="center", va="center", fontsize=7, color="gray")
        ax_exon.set_xlim(x[0], x[-1])
        ax_exon.set_ylim(0, 1)
        ax_exon.set_yticks([])
        ax_exon.set_xlabel("Position relative to this haplotype's own TSS (bp)")

    fig.suptitle(
        f"{result['gene']} -- {result['sample_id']} -- edited: {result['location_label']} -- "
        f"delta(log-odds {ctx.class_names[ctx.strong_idx]}/{ctx.class_names[ctx.weak_idx]}) = "
        f"{result['delta_log_odds']:+.4f}"
    )
    fig.tight_layout()
    out_path = ctx.fig_dir / f"knockout_{result['gene']}_{result['method']}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return fig


def plot_bulk_flip_rate_bars(gene_method_stats, bulk_genes, bulk_methods, method_colors, method_labels,
                              n_individuals, out_path):
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(bulk_genes))
    width = 0.8 / len(bulk_methods)
    for i, method in enumerate(bulk_methods):
        sub = gene_method_stats[gene_method_stats["method"] == method].set_index("gene").reindex(bulk_genes)
        offset = (i - (len(bulk_methods) - 1) / 2) * width
        ax.bar(x + offset, sub["flip_rate_pct"], width, label=method_labels[method], color=method_colors[method])
    ax.set_xticks(x)
    ax.set_xticklabels(bulk_genes, rotation=45, ha="right")
    ax.set_ylabel(f"% of {n_individuals} test individuals whose prediction flips")
    ax.set_title("Knockout-induced classification flip rate, by gene and location method")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.show()


def _logodds_to_prob(x):
    """sigmoid(log-odds) -> P(strong pigmentation). Exact inverse of log-odds = logit_strong -
    logit_weak, since for this binary softmax head that difference is the true log-odds."""
    return 1.0 / (1.0 + np.exp(-np.asarray(x, dtype=float)))


def _xaxis_view_labels(class_names, strong_idx, weak_idx):
    return {
        "log_odds": "log-odds",
        "probability": f"P({class_names[strong_idx]})",
    }


def plot_bulk_arrow_population_mean(gene_method_stats, bulk_methods, method_colors, method_labels,
                                     class_names, strong_idx, weak_idx, n_individuals, out_path):
    """Population-mean knockout arrow plot, one figure per x-axis mode ("log_odds" symlog, or
    "probability" of strong pigmentation in [0, 1]). Saved to `out_path` (log-odds) and a
    "_probability"-suffixed sibling, and returned as (gene_order_by_effect, figs) where `figs` is
    {"log_odds": fig, "probability": fig} for `render_carousel_with_toggle`."""
    gene_max_flip_rate = gene_method_stats.groupby("gene")["flip_rate"].max().sort_values(ascending=False)
    gene_order_by_effect = gene_max_flip_rate.index.tolist()
    n_methods = len(bulk_methods)
    row_height = 0.62

    figs = {}
    for x_mode in ("log_odds", "probability"):
        to_x = (lambda v: v) if x_mode == "log_odds" else _logodds_to_prob
        fig, ax = plt.subplots(figsize=(10, 0.62 * len(gene_order_by_effect) + 1.5))
        all_x = []
        for gi, gene in enumerate(gene_order_by_effect):
            y_base = (len(gene_order_by_effect) - 1 - gi) * row_height * n_methods
            for mi, method in enumerate(bulk_methods):
                row = gene_method_stats[(gene_method_stats["gene"] == gene) & (gene_method_stats["method"] == method)]
                if row.empty:
                    continue
                row = row.iloc[0]
                y = y_base + mi * row_height
                color = method_colors[method]
                x_baseline, x_perturbed = to_x(row["mean_baseline"]), to_x(row["mean_perturbed"])
                ax.annotate(
                    "", xy=(x_perturbed, y), xytext=(x_baseline, y),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=2.2, alpha=0.9, mutation_scale=14),
                )
                ax.scatter([x_baseline], [y], color=color, s=18, zorder=3, edgecolor="white", linewidth=0.4)
                all_x.extend([x_baseline, x_perturbed])
        ax.axvline(0 if x_mode == "log_odds" else 0.5, color="black", linewidth=0.7, linestyle="--", alpha=0.5)
        yticks = [(len(gene_order_by_effect) - 1 - gi) * row_height * n_methods + row_height * (n_methods - 1) / 2
                  for gi in range(len(gene_order_by_effect))]
        ax.set_yticks(yticks)
        ax.set_yticklabels(gene_order_by_effect)
        x_lo, x_hi = min(all_x), max(all_x)
        x_pad = 0.15 * (x_hi - x_lo)
        if x_mode == "log_odds":
            ax.set_xscale("symlog", linthresh=0.05)
            ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
            ax.set_xlabel(f"mean {class_names[strong_idx]}/{class_names[weak_idx]} log-odds, symlog scale  (dot = baseline, arrowhead = after knockout)")
        else:
            ax.set_xlim(max(0.0, x_lo - x_pad), min(1.0, x_hi + x_pad))
            ax.set_xlabel(f"mean P({class_names[strong_idx]})  (dot = baseline, arrowhead = after knockout)")
        ax.set_title(f"Population-mean knockout effect by gene x method (n={n_individuals} individuals)\ngenes ranked by max flip rate across methods")
        handles = [plt.Line2D([], [], color=method_colors[m], lw=2.2, label=method_labels[m]) for m in bulk_methods]
        ax.legend(handles=handles, loc="lower right", fontsize=9)
        ax.set_ylim(-row_height, len(gene_order_by_effect) * row_height * n_methods)
        fig.tight_layout()
        suffix = "" if x_mode == "log_odds" else "_probability"
        mode_out_path = out_path.parent / f"{out_path.stem}{suffix}{out_path.suffix}"
        fig.savefig(mode_out_path, dpi=140)
        plt.close(fig)
        figs[x_mode] = fig
    return gene_order_by_effect, figs


def plot_bulk_arrows_individual(bulk_df, genes, bulk_methods, method_labels, class_names,
                                 target_idx, other_idx, class_colors, strong_idx, weak_idx, out_dir):
    """Per-gene per-individual knockout arrows (baseline -> perturbed), one figure per gene per
    x-axis mode ("log_odds" or "probability" of strong pigmentation in [0, 1]), saved to `out_dir`
    and returned as (figures_by_view, labels) -- figures_by_view = {"log_odds": [...],
    "probability": [...]} -- for `render_carousel_with_toggle`. Individuals within each gene are
    always ordered by raw baseline log-odds (a monotonic transform doesn't change that order)."""
    true_label_colors = {
        class_names[target_idx]: class_colors[target_idx],
        class_names[other_idx]: class_colors[other_idx],
    }

    figures_by_view = {"log_odds": [], "probability": []}
    labels = []
    for gene in genes:
        gene_df = bulk_df[bulk_df["gene"] == gene]
        for x_mode in ("log_odds", "probability"):
            to_x = (lambda v: v) if x_mode == "log_odds" else _logodds_to_prob
            fig, axes = plt.subplots(1, len(bulk_methods), figsize=(15, 9), sharey=True)
            gene_x_lo = to_x(gene_df[["baseline_log_odds", "perturbed_log_odds"]].min().min())
            gene_x_hi = to_x(gene_df[["baseline_log_odds", "perturbed_log_odds"]].max().max())
            gene_x_pad = 0.08 * (gene_x_hi - gene_x_lo + 1e-6)
            for ax, method in zip(axes, bulk_methods):
                sub = gene_df[gene_df["method"] == method].sort_values("baseline_log_odds").reset_index(drop=True)
                for row_i, row in sub.iterrows():
                    color = true_label_colors.get(row["true_label"], "gray")
                    alpha = 0.9 if row["flipped"] else 0.25
                    lw = 1.6 if row["flipped"] else 0.7
                    x_baseline, x_perturbed = to_x(row["baseline_log_odds"]), to_x(row["perturbed_log_odds"])
                    ax.annotate(
                        "", xy=(x_perturbed, row_i), xytext=(x_baseline, row_i),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, alpha=alpha, mutation_scale=6),
                    )
                ax.set_xlim(gene_x_lo - gene_x_pad, gene_x_hi + gene_x_pad)
                ax.set_ylim(-1, len(sub))
                ax.axvline(0 if x_mode == "log_odds" else 0.5, color="black", linewidth=0.6, linestyle="--", alpha=0.4)
                n_flipped = int(sub["flipped"].sum())
                ax.set_title(f"{method_labels[method]}\n{n_flipped}/{len(sub)} flipped ({100 * n_flipped / len(sub):.0f}%)", fontsize=10)
                if x_mode == "log_odds":
                    ax.set_xlabel(f"{class_names[strong_idx]}/{class_names[weak_idx]} log-odds")
                else:
                    ax.set_xlabel(f"P({class_names[strong_idx]})")
            axes[0].set_ylabel("individuals, sorted by baseline log-odds")
            for ax in axes:
                ax.set_yticks([])
            handles = [plt.Line2D([], [], color=c, lw=1.6, label=lbl) for lbl, c in true_label_colors.items()]
            handles.append(plt.Line2D([], [], color="gray", lw=0.7, alpha=0.25, label="not flipped (faint)"))
            fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02), fontsize=9)
            fig.suptitle(f"{gene}: per-individual knockout arrows, baseline -> perturbed (n={sub['sample_id'].nunique()})", y=1.06)
            fig.tight_layout()
            suffix = "" if x_mode == "log_odds" else "_probability"
            out_path = out_dir / f"bulk_arrows_individual_{gene}{suffix}.png"
            fig.savefig(out_path, dpi=140, bbox_inches="tight")
            plt.close(fig)
            figures_by_view[x_mode].append(fig)
            print(f"saved {out_path}")
        labels.append(gene)

    return figures_by_view, labels
