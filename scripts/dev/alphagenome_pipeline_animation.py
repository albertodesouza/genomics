from pathlib import Path

import numpy as np
from manim import *


ASSET_DIR = Path(__file__).resolve().parents[2] / "assets" / "pipeline_animation" / "ppt_media"

BG = "#F7EFE3"
WHITE = "#3A2D26"
CYAN = "#7FA99B"
BLUE = "#6F8FA8"
PURPLE = "#A48A7B"
PINK = "#C76D57"
GREEN = "#8C9A64"
YELLOW = "#D6A34A"
ORANGE = "#C9824A"
MUTED = "#8F7D6D"


class AlphaGenomePipelineAnimation(Scene):
    """Conference animation for the AlphaGenome-derived superpopulation pipeline."""

    def construct(self):
        self.camera.background_color = BG
        self.add(self._background())
        self.title_scene()
        self.cohort_matrix_scene()
        self.alphagenome_scene()
        self.training_matrix_scene()
        self.normalization_scene()
        self.dataset_scene()
        self.closing_scene()

    def _background(self):
        group = VGroup()
        for radius, opacity, color in [(4.7, 0.10, BLUE), (3.3, 0.07, PURPLE), (2.0, 0.05, CYAN)]:
            group.add(Circle(radius=radius, color=color, stroke_width=2, stroke_opacity=opacity))
        for x in np.linspace(-6.8, 6.8, 11):
            group.add(Line([x, -4, 0], [x + 1.8, 4, 0], color=BLUE, stroke_opacity=0.055, stroke_width=1))
        for y in np.linspace(-3.6, 3.6, 7):
            group.add(Line([-7, y, 0], [7, y + 0.25, 0], color=PURPLE, stroke_opacity=0.06, stroke_width=1))
        return group

    def _maybe_image(self, name, height=2.0, opacity=0.22):
        path = ASSET_DIR / name
        if not path.exists():
            return VGroup()
        image = ImageMobject(str(path)).set_height(height).set_opacity(opacity)
        return image

    def _label(self, text, font_size=28, color=WHITE, weight=NORMAL):
        return Text(text, font_size=font_size, color=color, weight=weight)

    def _badge(self, text, color=CYAN, width=None):
        label = Text(text, font_size=22, color=WHITE, weight=BOLD)
        box_width = max(width or 0, label.width + 0.42)
        box_height = max(0.44, label.height + 0.18)
        box = RoundedRectangle(
            corner_radius=0.12,
            width=box_width,
            height=box_height,
            stroke_color=color,
            stroke_width=1.4,
            fill_color=color,
            fill_opacity=0.18,
        )
        label.move_to(box)
        return VGroup(box, label)

    def _flow_arrow(self, start, end, color=CYAN):
        return Arrow(start, end, buff=0.18, color=color, stroke_width=4, max_tip_length_to_length_ratio=0.12)

    def _signal_track(self, width=4.6, height=0.45, color=CYAN, phase=0.0, samples=90):
        xs = np.linspace(-width / 2, width / 2, samples)
        points = []
        for x in xs:
            y = 0.15 * np.sin(4.2 * x + phase) + 0.08 * np.sin(11.0 * x + phase)
            y += 0.18 * np.exp(-((x - 0.9) ** 2) / 0.18) + 0.12 * np.exp(-((x + 1.6) ** 2) / 0.35)
            points.append([x, np.clip(y, -height / 2, height / 2), 0])
        graph = VMobject(color=color, stroke_width=3)
        graph.set_points_smoothly(points)
        axis = Line([-width / 2, -height / 2, 0], [width / 2, -height / 2, 0], color=MUTED, stroke_opacity=0.4, stroke_width=1)
        return VGroup(axis, graph)

    def _strand_from_points(self, points, letters, color=CYAN, letter_offset=UP * 0.18, font_size=13, stroke_width=3):
        curve = VMobject(color=color, stroke_width=stroke_width)
        curve.set_points_smoothly(points)
        labels = VGroup()
        for i, (point, letter) in enumerate(zip(points, letters)):
            label_color = color if i % 5 else PINK
            offset = letter_offset(point, i) if callable(letter_offset) else letter_offset
            labels.add(Text(letter, font_size=font_size, color=label_color).move_to(np.array(point) + offset))
        strand = VGroup(curve, labels)
        strand.curve = curve
        strand.bases = labels
        return strand

    def _linear_strand(self, sequence, start, end, color=CYAN, letter_offset=RIGHT * 0.22, font_size=17, stroke_width=5):
        points = [interpolate(np.array(start), np.array(end), alpha) for alpha in np.linspace(0, 1, len(sequence))]
        return self._strand_from_points(points, sequence, color=color, letter_offset=letter_offset, font_size=font_size, stroke_width=stroke_width)

    def _dna(self, width=5.2, bases=18):
        pair_groups = VGroup()
        xs = np.linspace(-width / 2, width / 2, bases)
        pairs = [("A", "T"), ("C", "G"), ("G", "C"), ("T", "A")]
        top_points = []
        bot_points = []
        top_letters = []
        bottom_letters = []
        for i, x in enumerate(xs):
            y = 0.28 * np.sin(i * 0.9)
            top_points.append([x, y, 0])
            bot_points.append([x, -y, 0])
            top_base, bottom_base = pairs[i % len(pairs)]
            top_letters.append(top_base)
            bottom_letters.append(bottom_base)
            pair_groups.add(Line([x, y, 0], [x, -y, 0], color=BLUE, stroke_width=1.5, stroke_opacity=0.65))
        def top_outer_offset(point, _):
            return (UP if point[1] >= 0 else DOWN) * 0.16

        def bottom_outer_offset(point, _):
            return (DOWN if point[1] <= 0 else UP) * 0.16

        top = self._strand_from_points(top_points, top_letters, color=CYAN, letter_offset=top_outer_offset)
        bottom = self._strand_from_points(bot_points, bottom_letters, color=PURPLE, letter_offset=bottom_outer_offset)
        group = VGroup(top, bottom, pair_groups)
        group.top_strand = top
        group.bottom_strand = bottom
        group.strands = VGroup(top, bottom)
        group.base_pairs = pair_groups
        return group

    def _dna_growth(self, dna, run_time=1.5):
        base_reveals = []
        for i in range(len(dna.base_pairs)):
            base_reveals.append(FadeIn(dna.top_strand.bases[i], scale=0.9))
            base_reveals.append(FadeIn(dna.bottom_strand.bases[i], scale=0.9))
            base_reveals.append(FadeIn(dna.base_pairs[i], scale=0.9))
        return AnimationGroup(
            Create(VGroup(dna.top_strand.curve, dna.bottom_strand.curve), lag_ratio=0.0),
            LaggedStart(*base_reveals, lag_ratio=0.035),
            lag_ratio=0.0,
            run_time=run_time,
        )

    def _single_strand(self, width=4.6, bases=18, color=CYAN):
        group = VGroup()
        xs = np.linspace(-width / 2, width / 2, bases)
        letters = "ACGT"
        points = []
        for i, x in enumerate(xs):
            y = 0.16 * np.sin(i * 0.8)
            points.append([x, y, 0])
            base = Text(letters[i % 4], font_size=18, color=color if i % 5 else PINK).move_to([x, y + 0.22, 0])
            group.add(base)
        strand = VMobject(color=color, stroke_width=4)
        strand.set_points_smoothly(points)
        group.add_to_back(strand)
        group.strand = strand
        group.bases = VGroup(*group[1:])
        return group

    def _single_strand_growth(self, strand, run_time=1.2):
        return AnimationGroup(
            Create(strand.strand),
            LaggedStart(*[FadeIn(base, shift=RIGHT * 0.05) for base in strand.bases], lag_ratio=0.08),
            lag_ratio=0.0,
            run_time=run_time,
        )

    def _heatmap(self, rows=9, cols=18, cell=0.13):
        colors = [BLUE, CYAN, GREEN, YELLOW, PINK, PURPLE]
        grid = VGroup()
        for r in range(rows):
            for c in range(cols):
                val = (np.sin(r * 1.7 + c * 0.55) + 1) / 2
                color = colors[(r + c) % len(colors)]
                sq = Square(side_length=cell, stroke_width=0, fill_color=color, fill_opacity=0.20 + 0.65 * val)
                sq.move_to([(c - cols / 2) * cell, (rows / 2 - r) * cell, 0])
                grid.add(sq)
        frame = RoundedRectangle(
            corner_radius=0.06,
            width=cols * cell + 0.12,
            height=rows * cell + 0.12,
            stroke_color=CYAN,
            stroke_width=1.0,
            stroke_opacity=0.7,
        )
        frame.move_to(grid)
        return VGroup(grid, frame)

    def _heatmap_tile(self, rows=6, cols=14, cell=0.045, phase=0.0):
        grid = VGroup()
        for r in range(rows):
            for c in range(cols):
                val = (np.sin(r * 1.4 + c * 0.75 + phase) + 1) / 2
                color = interpolate_color(ManimColor(BLUE), ManimColor(PINK), val)
                square = Square(
                    side_length=cell,
                    stroke_width=0,
                    fill_color=color,
                    fill_opacity=0.35 + 0.60 * val,
                )
                square.move_to([(c - cols / 2) * cell, (rows / 2 - r) * cell, 0])
                grid.add(square)
        return grid

    def title_scene(self):
        title = Text("AlphaGenome-Derived Functional Signals", font_size=42, color=WHITE, weight=BOLD)
        subtitle = Text("Inferring superpopulation labels from individual DNA variation", font_size=24, color=MUTED)
        dna = self._dna(width=7.0, bases=28).scale(0.9).shift(DOWN * 1.0)
        group = VGroup(title, subtitle).arrange(DOWN, buff=0.22).shift(UP * 1.35)
        self.play(FadeIn(group, shift=UP * 0.25), self._dna_growth(dna), run_time=2.8)
        self.wait(3.6)
        self.play(FadeOut(group), FadeOut(dna), run_time=1.0)

    def cohort_matrix_scene(self):
        heading = self._label("Build cohort DNA windows before modeling", 32, WHITE, BOLD).to_edge(UP, buff=0.42)
        genes = ["DDB1", "EDAR", "HERC2", "MC1R", "MFSD12", "OCA2", "SLC24A5", "SLC45A2", "TCHH", "TYR", "TYRP1"]
        colors = [CYAN, BLUE, PURPLE, PINK, GREEN, YELLOW, ORANGE]
        chromosome_rows = [
            ("chr1", [("TCHH", 0.611)]),
            ("chr2", [("EDAR", 0.450)]),
            ("chr5", [("SLC45A2", 0.187)]),
            ("chr9", [("TYRP1", 0.092)]),
            ("chr11", [("DDB1", 0.454), ("TYR", 0.661)]),
            ("chr15", [("OCA2", 0.274), ("HERC2", 0.277), ("SLC24A5", 0.472)]),
            ("chr16", [("MC1R", 0.965)]),
            ("chr19", [("MFSD12", 0.061)]),
        ]
        gene_to_color = {gene: colors[i % len(colors)] for i, gene in enumerate(genes)}
        label_offsets = {
            "OCA2": LEFT * 0.30,
            "HERC2": RIGHT * 0.34 + UP * 0.04,
            "MC1R": LEFT * 0.30,
            "MFSD12": RIGHT * 0.34,
            "TYR": LEFT * 0.18,
        }
        window_offsets = {"MC1R": LEFT * 0.10}

        ref_label = Text("GRCh38 reference chromosomes", font_size=22, color=MUTED).move_to([-0.98, 2.12, 0])
        grch38 = VGroup()
        gene_marks = VGroup()
        gene_labels = VGroup()
        gene_ticks = VGroup()
        windows = VGroup()
        windows_by_gene = {}
        gene_to_chromosome = {}
        chrom_groups = {}
        for row, (chromosome, loci) in enumerate(chromosome_rows):
            y = 1.34 - row * 0.43
            label = Text(chromosome, font_size=15, color=MUTED, weight=BOLD).move_to([-4.55, y, 0])
            track = Line([-4.05, y, 0], [2.10, y, 0], color=MUTED, stroke_width=4.4, stroke_opacity=0.42)
            chrom_group = VGroup(label, track)
            for gene, position in loci:
                x = -4.05 + position * 6.15
                color = gene_to_color[gene]
                tick = Line([x, y - 0.10, 0], [x, y + 0.22, 0], color=color, stroke_width=3.4)
                gene_label = Text(gene, font_size=8, color=color).next_to(tick, UP, buff=0.05).shift(label_offsets.get(gene, ORIGIN))
                gene_ticks.add(tick)
                gene_labels.add(gene_label)
                gene_marks.add(VGroup(tick, gene_label))
                window = RoundedRectangle(
                    corner_radius=0.04,
                    width=0.40,
                    height=0.22,
                    stroke_color=YELLOW,
                    stroke_width=2,
                    fill_color=YELLOW,
                    fill_opacity=0.12,
                ).move_to([x, y, 0]).shift(window_offsets.get(gene, ORIGIN))
                windows.add(window)
                windows_by_gene[gene] = window
                gene_to_chromosome[gene] = chromosome
                chrom_group.add(tick, gene_label, window)
            chrom_groups[chromosome] = chrom_group
            grch38.add(chrom_group)
        size_bracket = VGroup(
            Line(LEFT * 0.20, RIGHT * 0.20, color=YELLOW, stroke_width=2.6),
            Line(LEFT * 0.20, LEFT * 0.20 + UP * 0.13, color=YELLOW, stroke_width=2.6),
            Line(RIGHT * 0.20, RIGHT * 0.20 + UP * 0.13, color=YELLOW, stroke_width=2.6),
            Text("512 kbp", font_size=13, color=YELLOW),
        ).arrange(DOWN, buff=0.03).next_to(windows_by_gene["TCHH"], DOWN, buff=0.08)
        left_reference = VGroup(ref_label, grch38, size_bracket)

        collection_title = VGroup(
            Text("1000 Genomes Project", font_size=21, color=WHITE, weight=BOLD),
            Text("high-coverage collection", font_size=16, color=MUTED),
        ).arrange(DOWN, buff=0.06).move_to([0.62, 2.05, 0])
        collection_box = RoundedRectangle(
            corner_radius=0.14,
            width=1.95,
            height=2.72,
            stroke_color=BLUE,
            stroke_width=1.4,
            fill_color=BLUE,
            fill_opacity=0.07,
        ).move_to([0.62, 0.10, 0])
        vcf_by_chromosome = {}
        vcf_cards = VGroup()
        vcf_by_chromosome = {}
        for chromosome, _ in chromosome_rows:
            rect = RoundedRectangle(corner_radius=0.08, width=1.34, height=0.23, stroke_color=BLUE, fill_color=BLUE, fill_opacity=0.14)
            txt = Text(f"{chromosome}.vcf.gz", font_size=10, color=WHITE)
            txt.move_to(rect)
            card = VGroup(rect, txt)
            vcf_cards.add(card)
            vcf_by_chromosome[chromosome] = card

        h1_title = Text("H1 windows by individual", font_size=18, color=PINK, weight=BOLD)
        output_frame = RoundedRectangle(
            corner_radius=0.12,
            width=3.55,
            height=3.00,
            stroke_color=CYAN,
            stroke_width=1.2,
            fill_color=CYAN,
            fill_opacity=0.045,
        ).move_to([4.42, 0.08, 0])
        individual_labels = VGroup()
        individual_xs = [3.05, 3.95, 4.85, 5.75]
        label_y = output_frame.get_top()[1] + 0.24
        for label, x in zip(["I1", "I2", "...", "I3202"], individual_xs):
            individual_labels.add(Text(label, font_size=16, color=YELLOW if label != "..." else MUTED, weight=BOLD).move_to([x, label_y, 0]))
        h1_title.next_to(output_frame, DOWN, buff=0.12)
        h1_rows = VGroup()
        h1_by_gene = {}
        gene_row_ys = np.linspace(1.35, -1.05, len(genes))
        for row, (gene, y) in enumerate(zip(genes, gene_row_ys)):
            row_group = VGroup()
            for col, x in enumerate(individual_xs):
                if col == 2:
                    row_group.add(Text("...", font_size=18, color=MUTED).move_to([x, y, 0]))
                    continue
                helix = self._dna(width=1.18, bases=6).scale(0.18)
                color_index = genes.index(gene)
                helix.top_strand.curve.set_color(colors[color_index % len(colors)])
                helix.bottom_strand.curve.set_color(colors[(color_index + 2) % len(colors)])
                helix.move_to([x, y, 0])
                row_group.add(helix)
            h1_by_gene[gene] = row_group
            h1_rows.add(row_group)
        h1_group = VGroup(output_frame, h1_title, individual_labels, h1_rows)

        split_colors = {"train": GREEN, "val": YELLOW, "test": PINK}
        split_circles = VGroup()
        split_targets = {
            "train": np.array([-4.05, -0.40, 0]),
            "val": np.array([-1.35, -0.40, 0]),
            "test": np.array([1.35, -0.40, 0]),
        }
        for split, center in split_targets.items():
            circle = Circle(radius=0.72, stroke_color=split_colors[split], stroke_width=2.0, fill_color=split_colors[split], fill_opacity=0.10).move_to(center)
            label = Text(split, font_size=21, color=split_colors[split], weight=BOLD).next_to(circle, UP, buff=0.12)
            split_circles.add(VGroup(circle, label))
        family_box = self._badge("family members grouped in the same split", ORANGE, width=4.85).move_to([0, -2.18, 0])
        def h1_individual_column(label, x, source_col):
            column = VGroup(individual_labels[source_col].copy())
            for row_group in h1_rows:
                column.add(row_group[source_col].copy())
            return column

        split_individual_columns = VGroup(
            h1_individual_column("I1", individual_xs[0], 0),
            h1_individual_column("I2", individual_xs[1], 1),
            h1_individual_column("I3202", individual_xs[3], 3),
        )
        split_assignment = ["train", "val", "test"]
        transit_split_columns = VGroup()
        for i, column in enumerate(split_individual_columns):
            split = split_assignment[i]
            source = column.get_center()
            for start_offset, end_offset in [
                (LEFT * 0.08 + DOWN * 0.04, LEFT * 0.22 + DOWN * 0.14),
                (RIGHT * 0.08 + DOWN * 0.04, RIGHT * 0.22 + DOWN * 0.14),
            ]:
                companion = column.copy().scale(0.82).move_to(source + start_offset).set_opacity(0.0)
                companion.target_split = split
                companion.target_offset = end_offset
                transit_split_columns.add(companion)

        matrix = VGroup()
        shown_genes = ["DDB1", "EDAR", "HERC2", "...", "TYRP1"]
        x_positions = [-1.55, -0.65, 0.25, 1.15, 2.05]
        matrix.add(Text("11 gene windows", font_size=19, color=YELLOW, weight=BOLD).move_to([0.55, 1.35, 0]))
        for i, gene in enumerate(shown_genes):
            matrix.add(Text(gene, font_size=14, color=colors[i % len(colors)] if gene != "..." else MUTED).move_to([x_positions[i], 0.98, 0]))
        split_rows = [("train", 0.45, GREEN, "Many rows"), ("val", -0.35, YELLOW, "Held-Out"), ("test", -1.15, PINK, "Held-Out")]
        for row, (split, y, color, note) in enumerate(split_rows):
            band = RoundedRectangle(corner_radius=0.10, width=6.15, height=0.58, stroke_color=color, stroke_width=1.0, fill_color=color, fill_opacity=0.075).move_to([0.0, y, 0])
            matrix.add(band)
            matrix.add(Text(split, font_size=17, color=color, weight=BOLD).move_to([-2.65, y, 0]))
            matrix.add(Text(note, font_size=11, color=MUTED).move_to([-2.65, y - 0.20, 0]))
            for col, gene in enumerate(shown_genes):
                if gene == "...":
                    matrix.add(Text("...", font_size=22, color=MUTED).move_to([x_positions[col], y, 0]))
                    continue
                cell = self._dna(width=1.15, bases=5).scale(0.17)
                cell.top_strand.curve.set_color(colors[col % len(colors)])
                cell.bottom_strand.curve.set_color(colors[(col + row + 2) % len(colors)])
                cell.move_to([x_positions[col], y, 0])
                matrix.add(cell)
        matrix_frame = RoundedRectangle(corner_radius=0.12, width=6.55, height=3.08, stroke_color=CYAN, stroke_width=1.3, fill_color=CYAN, fill_opacity=0.05).move_to([0.0, -0.02, 0])
        matrix.add_to_back(matrix_frame)
        matrix_note = self._badge("split × individual × gene × double-stranded H1 window", CYAN, width=6.25).next_to(matrix, DOWN, buff=0.18)
        matrix_count = Text("3,202 individuals; family members never cross splits", font_size=19, color=MUTED).next_to(matrix_note, DOWN, buff=0.14)

        self.play(Write(heading), FadeIn(ref_label, shift=UP * 0.12), LaggedStart(*[Create(chrom_groups[chromosome][1]) for chromosome, _ in chromosome_rows], lag_ratio=0.05), run_time=1.4)
        self.play(LaggedStart(*[FadeIn(chrom_groups[chromosome][0], shift=RIGHT * 0.06) for chromosome, _ in chromosome_rows], lag_ratio=0.04), run_time=0.6)
        self.play(LaggedStart(*[FadeIn(mark, shift=UP * 0.08) for mark in gene_marks], lag_ratio=0.05), run_time=1.5)
        self.play(LaggedStart(*[Create(window) for window in windows], lag_ratio=0.04), FadeIn(size_bracket, shift=DOWN * 0.06), gene_ticks.animate.set_opacity(0.70), run_time=1.1)
        self.wait(2.0)
        self.play(
            left_reference.animate.scale(0.74).move_to([-4.18, 0.00, 0]),
            run_time=1.1,
        )
        for row, (chromosome, _) in enumerate(chromosome_rows):
            vcf_by_chromosome[chromosome].move_to([0.62, 1.18 - row * 0.31, 0])
        self.play(
            FadeIn(collection_box),
            FadeIn(collection_title, shift=LEFT * 0.10),
            LaggedStart(*[FadeIn(card, shift=RIGHT * 0.14) for card in vcf_cards], lag_ratio=0.07),
            run_time=1.2,
        )
        self.play(
            FadeIn(output_frame),
            FadeIn(h1_title, shift=UP * 0.05),
            FadeIn(individual_labels, shift=DOWN * 0.08),
            run_time=0.9,
        )
        for gene in genes:
            chromosome = gene_to_chromosome[gene]
            moving = VGroup(windows_by_gene[gene].copy())
            self.add(moving)
            self.play(
                windows_by_gene[gene].animate.set_opacity(0.20),
                moving.animate.move_to(collection_box.get_center()).scale(0.34).set_opacity(0.0),
                FadeIn(h1_by_gene[gene], shift=DOWN * 0.06),
                run_time=0.35,
            )
            self.remove(moving)
        self.wait(1.6)
        self.play(FadeOut(VGroup(left_reference, collection_box, collection_title, vcf_cards)), run_time=0.8)
        self.play(FadeIn(split_circles, shift=UP * 0.14), FadeIn(family_box, shift=UP * 0.10), run_time=1.0)
        self.add(split_individual_columns, transit_split_columns)
        visible_h1_rows = VGroup(*[h1_by_gene[gene] for gene in genes])
        self.remove(visible_h1_rows, individual_labels)
        self.play(
            *[
                column.animate.move_to(split_targets[split_assignment[i]] + UP * 0.10).scale(0.44).set_opacity(0.80)
                for i, column in enumerate(split_individual_columns)
            ],
            *[
                column.animate.move_to(split_targets[column.target_split] + column.target_offset).scale(0.40).set_opacity(0.58)
                for column in transit_split_columns
            ],
            run_time=1.6,
        )
        self.play(FadeOut(VGroup(output_frame, h1_title)), run_time=0.6)
        self.wait(2.2)
        remaining = VGroup(heading, split_circles, family_box, split_individual_columns, transit_split_columns)
        self.play(FadeOut(remaining), run_time=0.8)
        self.remove(
            left_reference,
            collection_box,
            collection_title,
            vcf_cards,
            output_frame,
            h1_title,
            h1_rows,
            individual_labels,
            split_circles,
            family_box,
            split_individual_columns,
            transit_split_columns,
            heading,
        )

    def alphagenome_scene(self):
        heading = self._label("AlphaGenome: one strand in, 5930 human tracks out", 32, WHITE, BOLD).to_edge(UP, buff=0.45)
        helix = self._dna(width=2.9, bases=12).scale(0.58).shift(LEFT * 5.45)
        plus_sequence = "ACGTACGTACGT"
        minus_sequence = "TGCATGCATGCA"
        plus = self._linear_strand(plus_sequence, UP * 0.95, DOWN * 0.95, color=CYAN, letter_offset=RIGHT * 0.14, font_size=8, stroke_width=3).shift(LEFT * 4.10 + UP * 1.25)
        minus = self._linear_strand(minus_sequence, UP * 0.95, DOWN * 0.95, color=PURPLE, letter_offset=RIGHT * 0.14, font_size=8, stroke_width=3).shift(LEFT * 4.10 + DOWN * 1.25)
        models = VGroup()
        for y in [1.25, -1.25]:
            model = RoundedRectangle(corner_radius=0.18, width=2.55, height=1.1, stroke_color=CYAN, fill_color=PURPLE, fill_opacity=0.18)
            model_text = VGroup(Text("AlphaGenome", font_size=25, color=WHITE, weight=BOLD)).arrange(DOWN, buff=0.10)
            model_text.move_to(model)
            model_group = VGroup(model, model_text).shift(LEFT * 1.70 + UP * y)
            models.add(model_group)
        top_tracks = VGroup()
        bottom_tracks = VGroup()
        track_colors = [CYAN, GREEN, YELLOW, PINK, PURPLE, ORANGE, BLUE]
        for i, color in enumerate(track_colors):
            top_tracks.add(self._signal_track(width=2.35, height=0.24, color=color, phase=i).scale(0.58).shift(RIGHT * 1.35 + UP * (2.05 - i * 0.19)))
            bottom_tracks.add(self._signal_track(width=2.35, height=0.24, color=color, phase=i + 0.7).scale(0.58).shift(RIGHT * 1.35 + DOWN * (0.45 + i * 0.19)))
        top_count = self._badge("5930 tracks", ORANGE, width=1.62).next_to(top_tracks, DOWN, buff=0.14)
        bottom_count = self._badge("5930 tracks", ORANGE, width=1.62).next_to(bottom_tracks, DOWN, buff=0.14)
        ellipses = VGroup(Text("...", font_size=26, color=MUTED).next_to(top_tracks[-1], DOWN, buff=0.01), Text("...", font_size=26, color=MUTED).next_to(bottom_tracks[-1], DOWN, buff=0.01))
        selected_indices = [1, 3, 5]
        selected_specs = [
            ("CL:1000458 RNA-Seq", GREEN),
            ("CL:0000346 RNA-Seq", PINK),
            ("CL:2000092 RNA-Seq", ORANGE),
        ]
        selected_labels = VGroup()
        selected_targets = []
        for strand, y_start in [("top", 1.48), ("bottom", -0.50)]:
            for i, (label, color) in enumerate(selected_specs):
                target = np.array([3.5, y_start - i * 0.43, 0])
                selected_targets.append((strand, selected_indices[i], target))
                text = Text(label, font_size=11, color=color, weight=BOLD).move_to([5.25, target[1], 0])
                selected_labels.add(text)
        selected_heatmap = self._heatmap_tile(rows=6, cols=22, cell=0.056, phase=1.1).move_to([4.72, 0.42, 0])
        heatmap_label = self._badge("6 selected tracks as one heatmap", GREEN, width=3.18).scale(0.64).next_to(selected_heatmap, DOWN, buff=0.18)
        heatmap_group = VGroup(selected_heatmap, heatmap_label)
        self.play(Write(heading), self._dna_growth(helix), FadeIn(models, scale=0.96), run_time=1.8)
        self.play(
            FadeOut(helix.base_pairs),
            helix.top_strand.animate.shift(UP * 0.24),
            helix.bottom_strand.animate.shift(DOWN * 0.24),
            run_time=0.8,
        )
        self.play(
            Transform(helix.top_strand, plus),
            Transform(helix.bottom_strand, minus),
            run_time=1.6,
        )
        plus = helix.top_strand
        minus = helix.bottom_strand
        self.play(plus.animate.move_to(models[0]).set_opacity(0.0), minus.animate.move_to(models[1]).set_opacity(0.0), run_time=1.3)
        self.play(
            LaggedStart(*[Create(t) for t in top_tracks], lag_ratio=0.08),
            LaggedStart(*[Create(t) for t in bottom_tracks], lag_ratio=0.08),
            FadeIn(ellipses),
            FadeIn(top_count, shift=UP * 0.12),
            FadeIn(bottom_count, shift=UP * 0.12),
            run_time=2.2,
        )
        self.wait(1.4)
        self.play(
            *[
                track.animate.set_color(MUTED).set_opacity(0.16)
                for i, track in enumerate(top_tracks)
                if i not in selected_indices
            ],
            *[
                track.animate.set_color(MUTED).set_opacity(0.16)
                for i, track in enumerate(bottom_tracks)
                if i not in selected_indices
            ],
            ellipses.animate.set_opacity(0.24),
            top_count.animate.set_opacity(0.35),
            bottom_count.animate.set_opacity(0.35),
            FadeIn(selected_labels, shift=LEFT * 0.12),
            run_time=1.5,
        )
        self.play(
            *[
                (top_tracks[index] if strand == "top" else bottom_tracks[index]).animate.move_to(target).scale(1.12)
                for strand, index, target in selected_targets
            ],
            run_time=1.5,
        )
        self.wait(2.2)
        self.play(
            ReplacementTransform(
                VGroup(
                    *[top_tracks[index].copy() for index in selected_indices],
                    *[bottom_tracks[index].copy() for index in selected_indices],
                ),
                selected_heatmap,
            ),
            FadeIn(heatmap_label, shift=UP * 0.06),
            FadeOut(VGroup(top_tracks, bottom_tracks, selected_labels, ellipses, top_count, bottom_count)),
            run_time=1.6,
        )
        self.wait(2.6)
        self.play(FadeOut(VGroup(models, heatmap_group, heading)), run_time=0.8)

    def track_selection_scene(self):
        heading = self._label("Filter 5930 human tracks to three RNA-seq tissue ontologies", 30, WHITE, BOLD).to_edge(UP, buff=0.45)
        selected = ["CL:1000458", "CL:0000346", "CL:2000092"]
        tracks = VGroup()
        colors = [MUTED, CYAN, MUTED, GREEN, MUTED, YELLOW, MUTED, MUTED, MUTED]
        for i, color in enumerate(colors):
            opacity = 1.0 if color != MUTED else 0.18
            tr = self._signal_track(width=5.0, height=0.36, color=color, phase=i).shift(UP * (1.85 - i * 0.38))
            tr.set_opacity(opacity)
            tracks.add(tr)
        source = self._badge("5930 human AlphaGenome tracks", ORANGE, width=3.6).shift(DOWN * 1.45)
        labels = VGroup(*[self._badge(s, [CYAN, GREEN, YELLOW][i], width=1.6) for i, s in enumerate(selected)]).arrange(RIGHT, buff=0.35).shift(DOWN * 2.08)
        self.play(Write(heading), LaggedStart(*[Create(t) for t in tracks], lag_ratio=0.06), run_time=1.5)
        self.play(FadeIn(source, shift=UP * 0.2), run_time=0.6)
        self.play(tracks.animate.scale(0.92), source.animate.set_opacity(0.22), FadeIn(labels, shift=UP * 0.2), run_time=1.3)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, tracks, source, labels)), run_time=0.7)

    def training_matrix_scene(self):
        heading = self._label("Build the training matrix from strand-specific AlphaGenome tracks", 29, WHITE, BOLD).to_edge(UP, buff=0.42)
        genes = ["DDB1", "EDAR", "HERC2", "MC1R", "MFSD12", "OCA2", "SLC24A5", "SLC45A2", "TCHH", "TYR", "TYRP1"]
        individuals = ["I1", "I2", "...", "I3202"]
        track_colors = [GREEN, PINK, ORANGE, GREEN, PINK, ORANGE]

        def tiny_track(width=0.44, height=0.035, color=CYAN, phase=0.0):
            xs = np.linspace(-width / 2, width / 2, 16)
            points = [[x, 0.45 * height * np.sin(12 * x + phase), 0] for x in xs]
            line = VMobject(color=color, stroke_width=1.05)
            line.set_points_smoothly(points)
            return line

        def six_track_stack(width=0.44, colors=track_colors, phase=0.0, vbuff=0.020):
            stack = VGroup()
            for i, color in enumerate(colors):
                stack.add(tiny_track(width=width, color=color, phase=phase + i * 0.55).shift(DOWN * i * vbuff))
            stack.arrange(DOWN, buff=0.010)
            return stack

        left_frame = RoundedRectangle(corner_radius=0.12, width=3.05, height=3.68, stroke_color=CYAN, stroke_width=1.2, fill_color=CYAN, fill_opacity=0.045).move_to([-4.75, -0.08, 0])
        left_title = self._badge("H1 double-stranded windows", CYAN, width=2.80).scale(0.70).next_to(left_frame, UP, buff=0.10)
        gene_labels = VGroup()
        input_columns = VGroup()
        col_xs = [-5.35, -4.85, -4.35, -3.85]
        row_ys = np.linspace(1.32, -1.65, len(genes))
        for gene, y in zip(genes, row_ys):
            gene_labels.add(Text(gene, font_size=9, color=MUTED).move_to([-6.02, y, 0]))
        for col, (individual, x) in enumerate(zip(individuals, col_xs)):
            column = VGroup(Text(individual, font_size=13, color=YELLOW if individual != "..." else MUTED, weight=BOLD).move_to([x, 1.6, 0]))
            for row, y in enumerate(row_ys):
                if individual == "...":
                    column.add(Text("...", font_size=10, color=MUTED).move_to([x, y, 0]))
                    continue
                helix = self._dna(width=0.95, bases=5).scale(0.135).move_to([x, y, 0])
                helix.top_strand.curve.set_color(track_colors[row % 3])
                helix.bottom_strand.curve.set_color(track_colors[(row + 1) % 3])
                column.add(helix)
            input_columns.add(column)
        input_group = VGroup(left_frame, left_title, gene_labels, input_columns)

        models = VGroup()
        for y, color, label in [(0.70, BLUE, "+ strand"), (-0.35, PURPLE, "- strand")]:
            block = RoundedRectangle(corner_radius=0.14, width=1.62, height=0.58, stroke_color=color, fill_color=color, fill_opacity=0.12).move_to([-1.05, y, 0])
            text = VGroup(Text("AlphaGenome", font_size=14, color=WHITE, weight=BOLD), Text(label, font_size=10, color=MUTED)).arrange(DOWN, buff=0.04).move_to(block)
            models.add(VGroup(block, text))

        output_frame = RoundedRectangle(corner_radius=0.12, width=3.62, height=3.72, stroke_color=GREEN, stroke_width=1.2, fill_color=GREEN, fill_opacity=0.04).move_to([3.95, -0.06, 0])
        output_title = self._badge("Training matrix: compact image columns", GREEN, width=3.55).scale(0.68).next_to(output_frame, UP, buff=0.10)
        output_columns = VGroup()
        out_xs = [2.78, 3.56, 4.34, 5.12]
        row_gap = 0.265
        for individual, x in zip(individuals, out_xs):
            column = VGroup(Text(individual, font_size=13, color=YELLOW if individual != "..." else MUTED, weight=BOLD).move_to([x, 1.6, 0]))
            for row in range(len(genes)):
                y = 1.10 - row * row_gap
                if individual == "...":
                    tile = RoundedRectangle(corner_radius=0.025, width=0.48, height=0.18, stroke_color=MUTED, stroke_width=0.7, fill_color=MUTED, fill_opacity=0.08).move_to([x, y, 0])
                    dots = Text("...", font_size=9, color=MUTED).move_to(tile)
                    image_tile = VGroup(tile, dots)
                else:
                    image_tile = self._heatmap_tile(rows=6, cols=12, cell=0.033, phase=row * 0.55 + x).move_to([x, y, 0])
                image_tile.set_opacity(0.0)
                column.add(image_tile)
            output_columns.add(column)
        track_note = self._badge("11 heatmap tiles: each tile encodes 6 selected tracks", ORANGE, width=4.25).scale(0.62).next_to(output_frame, DOWN, buff=0.10)
        output_group = VGroup(output_frame, output_title, output_columns, track_note)

        self.play(Write(heading), FadeIn(input_group, shift=RIGHT * 0.10), run_time=1.3)
        self.play(FadeIn(models, scale=0.96), FadeIn(VGroup(output_frame, output_title), shift=LEFT * 0.10), run_time=1.0)
        self.play(FadeIn(VGroup(*[column[0] for column in output_columns])), FadeIn(track_note, shift=UP * 0.08), run_time=0.7)
        self.wait(1.0)

        def output_tile_for(col, row):
            return output_columns[col][row + 1]

        def animate_helix(col, row, detailed=False):
            source = input_columns[col][row + 1].copy()
            source.generate_target()
            source.target.scale(2.5 if detailed else 1.55).move_to([-2.45, 0.18, 0])
            tile = output_tile_for(col, row)
            gene_label = gene_labels[row]
            focus_label = self._badge(f"{individuals[col]} / {genes[row]}", YELLOW, width=1.70).scale(0.62).move_to([-2.45, -1.26, 0])
            self.play(gene_label.animate.set_color(YELLOW).scale(1.18), MoveToTarget(source), FadeIn(focus_label, shift=UP * 0.05), run_time=0.8 if detailed else 0.22)
            self.play(
                FadeOut(source.base_pairs),
                source.top_strand.animate.shift(UP * (0.32 if detailed else 0.18)),
                source.bottom_strand.animate.shift(DOWN * (0.32 if detailed else 0.18)),
                run_time=0.6 if detailed else 0.14,
            )
            plus = self._linear_strand("ACGTAC", UP * 0.30, DOWN * 0.30, color=CYAN, font_size=7, stroke_width=2).scale(0.74 if detailed else 0.50).move_to(models[0])
            minus = self._linear_strand("TGCATG", UP * 0.30, DOWN * 0.30, color=PURPLE, font_size=7, stroke_width=2).scale(0.74 if detailed else 0.50).move_to(models[1])
            self.play(Transform(source.top_strand, plus), Transform(source.bottom_strand, minus), run_time=0.8 if detailed else 0.16)
            self.play(source.top_strand.animate.move_to(models[0]).set_opacity(0.0), source.bottom_strand.animate.move_to(models[1]).set_opacity(0.0), run_time=0.5 if detailed else 0.10)
            heatmap_burst = self._heatmap_tile(rows=6, cols=14, cell=0.055 if detailed else 0.037, phase=row + col).move_to([0.95, 0.18, 0])
            burst_label = Text("heatmap", font_size=12 if detailed else 8, color=GREEN).next_to(heatmap_burst, DOWN, buff=0.05)
            self.play(FadeIn(VGroup(heatmap_burst, burst_label), shift=RIGHT * 0.10), run_time=0.5 if detailed else 0.10)
            self.play(FadeOut(heatmap_burst, scale=0.75), tile.animate.set_opacity(1.0), FadeOut(VGroup(source, focus_label, burst_label)), run_time=0.8 if detailed else 0.18)
            self.play(gene_label.animate.set_color(MUTED).scale(1 / 1.18), run_time=0.2 if detailed else 0.07)

        animate_helix(0, 0, detailed=True)
        animate_helix(0, 1, detailed=True)

        remaining_steps = [(0, row) for row in range(2, len(genes))]
        remaining_steps += [(col, row) for col in [1, 3] for row in range(len(genes))]
        remaining_sources = VGroup()
        for col, row in remaining_steps:
            source = input_columns[col][row + 1].copy()
            source.scale(0.72).set_opacity(0.72)
            remaining_sources.add(source)
        remaining_sources.arrange_in_grid(rows=5, cols=7, buff=(0.07, 0.07)).move_to([-2.42, 0.12, 0])
        remaining_tiles = VGroup()
        for col, row in remaining_steps:
            remaining_tiles.add(output_tile_for(col, row))
        plus_cloud = VGroup(*[
            self._linear_strand("ACGTAC", LEFT * 0.18, RIGHT * 0.18, color=CYAN, font_size=4, stroke_width=1.1).scale(0.35)
            for _ in remaining_steps
        ]).arrange_in_grid(rows=5, cols=7, buff=(0.045, 0.030)).move_to(models[0])
        minus_cloud = VGroup(*[
            self._linear_strand("TGCATG", LEFT * 0.18, RIGHT * 0.18, color=PURPLE, font_size=4, stroke_width=1.1).scale(0.35)
            for _ in remaining_steps
        ]).arrange_in_grid(rows=5, cols=7, buff=(0.045, 0.030)).move_to(models[1])
        heatmap_cloud = VGroup(*[
            self._heatmap_tile(rows=6, cols=10, cell=0.021, phase=i * 0.33)
            for i, _ in enumerate(remaining_steps)
        ]).arrange_in_grid(rows=5, cols=7, buff=(0.045, 0.030)).move_to([0.98, 0.12, 0])
        self.play(
            LaggedStart(*[gene_labels[row].animate.set_color(YELLOW) for _, row in remaining_steps[:9]], lag_ratio=0.03),
            FadeIn(remaining_sources, shift=RIGHT * 0.12),
            run_time=1.0,
        )
        self.play(
            FadeOut(remaining_sources, scale=0.78),
            FadeIn(plus_cloud, shift=UP * 0.10),
            FadeIn(minus_cloud, shift=DOWN * 0.10),
            run_time=0.8,
        )
        self.play(
            FadeOut(VGroup(plus_cloud, minus_cloud), scale=0.70),
            FadeIn(heatmap_cloud, scale=0.95),
            run_time=0.8,
        )
        self.play(
            FadeOut(heatmap_cloud, scale=0.85),
            remaining_tiles.animate.set_opacity(1.0),
            VGroup(*output_columns[2][1:]).animate.set_opacity(0.35),
            run_time=1.2,
        )
        self.play(gene_labels.animate.set_color(MUTED), run_time=0.25)

        scene_group = VGroup(heading, input_group, models, output_group)
        self.wait(2.4)
        self.play(scene_group.animate.set_opacity(0.0), run_time=0.7)
        self.remove(scene_group, heatmap_cloud, remaining_sources, plus_cloud, minus_cloud)

    def normalization_scene(self):
        heading = self._label("Crop each heatmap to 32 kbp and log-normalize per track", 29, WHITE, BOLD).to_edge(UP, buff=0.45)
        def heatmap_stack(cols, cell, center, phases, x_label, color):
            stack = VGroup()
            for i, phase in enumerate(phases):
                # Visually compressed grid; labels preserve the true tensor shape.
                tile = self._heatmap_tile(rows=33, cols=cols, cell=cell, phase=phase)
                tile.set_opacity(0.92 - i * 0.16)
                tile.shift(RIGHT * i * 0.16 + UP * i * 0.10)
                stack.add(tile)
            stack.move_to(center)
            label = self._badge(x_label, color, width=3.05).scale(0.70).next_to(stack, UP, buff=0.18)
            return VGroup(stack, label)

        raw = heatmap_stack(32, 0.029, [-3.65, 0.30, 0], [0.2, 1.1, 2.0], "Individuals × 66 × 512k", CYAN)
        crop = Rectangle(width=0.86, height=1.28, stroke_color=YELLOW, stroke_width=3.0, fill_color=YELLOW, fill_opacity=0.08).move_to(raw[0])
        crop_label = Text("center 32 kbp", font_size=19, color=YELLOW, weight=BOLD).next_to(crop, DOWN, buff=0.15)
        cropped = heatmap_stack(16, 0.041, [0.00, 0.30, 0], [0.7, 1.6, 2.5], "Individuals × 66 × 32k", YELLOW)
        normalized = heatmap_stack(16, 0.041, [3.55, 0.30, 0], [1.2, 2.1, 3.0], "Same shape after normalization", GREEN)
        track_guides = VGroup()
        for i in range(6):
            y = normalized[0][0].get_center()[1] + 0.55 - i * 0.18
            track_guides.add(Line([3.00, y, 0], [4.24, y, 0], color=interpolate_color(ManimColor(BLUE), ManimColor(PINK), i / 5), stroke_width=1.6, stroke_opacity=0.65))
        formula = MathTex(r"x_{norm}=\frac{\ln(x+1)}{\ln(x_{max,track}+1)}", font_size=36, color=WHITE).shift(DOWN * 1.45)
        self.play(Write(heading), FadeIn(raw, shift=UP * 0.10), run_time=1.4)
        self.wait(0.8)
        self.play(Create(crop), FadeIn(crop_label), run_time=1.0)
        self.wait(1.0)
        self.play(FadeIn(cropped, shift=RIGHT * 0.15), raw.animate.set_opacity(0.28), run_time=1.2)
        self.wait(1.0)
        self.play(FadeIn(normalized, shift=RIGHT * 0.15), cropped.animate.set_opacity(0.55), FadeIn(track_guides), Write(formula), run_time=1.4)
        self.wait(6.0)
        self.play(FadeOut(VGroup(heading, raw, crop, crop_label, cropped, normalized, track_guides, formula)), run_time=0.8)

    def dataset_scene(self):
        heading = self._label("Train and evaluate a CNN", 32, WHITE, BOLD).to_edge(UP, buff=0.45)
        labels = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        label_colors = [PINK, ORANGE, YELLOW, CYAN, GREEN]

        heatmap_tiles = VGroup(*[self._heatmap_tile(rows=6, cols=14, cell=0.038, phase=0.4 + i * 0.48) for i in range(11)])
        heatmap_tiles.arrange(DOWN, buff=0.018)
        heatmap_frame = RoundedRectangle(
            corner_radius=0.12,
            width=heatmap_tiles.width + 0.24,
            height=heatmap_tiles.height + 0.24,
            stroke_color=CYAN,
            stroke_width=1.5,
            fill_color=CYAN,
            fill_opacity=0.045,
        ).move_to(heatmap_tiles)
        input_card = VGroup(heatmap_frame, heatmap_tiles).scale(1.02).move_to([-4.20, -0.10, 0])
        input_label = VGroup(
            Text("held-out individual", font_size=19, color=MUTED),
            Text("66 × 32768", font_size=23, color=CYAN, weight=BOLD),
        ).arrange(DOWN, buff=0.08).next_to(input_card, DOWN, buff=0.18)

        cnn_box = RoundedRectangle(corner_radius=0.22, width=2.10, height=1.35, stroke_color=PINK, stroke_width=2.0, fill_color=PINK, fill_opacity=0.13).move_to([-0.70, -0.05, 0])
        cnn_text = Text("CNN", font_size=38, color=WHITE, weight=BOLD).move_to(cnn_box)
        cnn_note = Text("Learned Filters → Class Scores", font_size=15, color=MUTED).next_to(cnn_box, DOWN, buff=0.16)
        cnn_group = VGroup(cnn_box, cnn_text, cnn_note)

        probabilities = VGroup()
        probability_values = [("AFR", 0.07, PINK), ("AMR", 0.10, ORANGE), ("EAS", 0.08, YELLOW), ("EUR", 0.61, CYAN), ("SAS", 0.14, GREEN)]
        for label, value, color in probability_values:
            name = Text(label, font_size=18, color=color, weight=BOLD)
            rail = RoundedRectangle(corner_radius=0.04, width=1.85, height=0.14, stroke_width=0, fill_color=MUTED, fill_opacity=0.16)
            bar = RoundedRectangle(corner_radius=0.04, width=1.85 * value / 0.61, height=0.14, stroke_width=0, fill_color=color, fill_opacity=0.86).align_to(rail, LEFT)
            value_text = Text(f"{value:.2f}", font_size=15, color=WHITE if label == "EUR" else MUTED)
            row = VGroup(name, VGroup(rail, bar), value_text).arrange(RIGHT, buff=0.14)
            probabilities.add(row)
        probabilities.arrange(DOWN, aligned_edge=LEFT, buff=0.16).move_to([3.05, 0.10, 0])
        prob_title = Text("Superpopulation Probabilities", font_size=20, color=WHITE, weight=BOLD).next_to(probabilities, UP, buff=0.22)
        comparison = VGroup(
            VGroup(Text("top", font_size=16, color=MUTED), Text("EUR", font_size=28, color=CYAN, weight=BOLD)).arrange(DOWN, buff=0.03),
            Text("=", font_size=28, color=MUTED),
            VGroup(Text("true", font_size=16, color=MUTED), Text("EUR", font_size=28, color=GREEN, weight=BOLD)).arrange(DOWN, buff=0.03),
        ).arrange(RIGHT, buff=0.28).move_to([3.05, -1.75, 0])
        match_box = RoundedRectangle(corner_radius=0.10, width=1.18, height=0.38, stroke_color=GREEN, stroke_width=1.1, fill_color=GREEN, fill_opacity=0.14).next_to(comparison, RIGHT, buff=0.28)
        match_text = Text("match", font_size=18, color=GREEN, weight=BOLD).move_to(match_box)
        match_group = VGroup(match_box, match_text)

        architecture_title = Text("CNN architecture", font_size=17, color=MUTED, weight=BOLD).move_to([0, -2.32, 0])
        architecture_steps = [
            ("Input\n66 × 32768", CYAN, 0.96),
            ("Conv1\n16 × 11 × 1024", PINK, 1.05),
            ("Conv2\n32 × 11 × 512", PINK, 1.05),
            ("Conv3\n64 × 11 × 256", PINK, 1.05),
            ("Pool\n64 × 11 × 1", BLUE, 0.92),
            ("Flatten\n704", PURPLE, 0.72),
            ("FC\n256", PURPLE, 0.62),
            ("Output\n5 classes", ORANGE, 0.82),
        ]
        architecture = VGroup()
        for text, color, width in architecture_steps:
            box = RoundedRectangle(corner_radius=0.06, width=width, height=0.42, stroke_color=color, stroke_width=0.8, fill_color=color, fill_opacity=0.075)
            label = Text(text, font_size=8.5, color=WHITE, weight=BOLD, line_spacing=0.72).move_to(box)
            architecture.add(VGroup(box, label))
        architecture.arrange(RIGHT, buff=0.08).move_to([0, -2.82, 0])
        architecture_lines = VGroup()
        for left, right in zip(architecture[:-1], architecture[1:]):
            architecture_lines.add(Line(left.get_right(), right.get_left(), color=MUTED, stroke_width=1.2, stroke_opacity=0.55))
        architecture_group = VGroup(architecture_title, architecture_lines, architecture)

        moving_heatmap = input_card.copy()
        self.play(Write(heading), run_time=1.2)
        self.play(FadeIn(input_card, shift=UP * 0.10), FadeIn(input_label, shift=UP * 0.06), FadeIn(cnn_group, scale=0.96), FadeIn(architecture_group, shift=UP * 0.08), run_time=1.4)
        self.wait(3.0)
        self.add(moving_heatmap)
        self.play(moving_heatmap.animate.scale(0.28).move_to(cnn_box).set_opacity(0.0), cnn_box.animate.set_fill(PINK, opacity=0.22), run_time=1.2)
        self.play(FadeIn(VGroup(prob_title, probabilities), shift=LEFT * 0.12), run_time=1.1)
        self.wait(2.0)
        self.play(FadeIn(comparison, shift=UP * 0.10), FadeIn(match_group, shift=UP * 0.10), run_time=0.8)
        self.wait(5.5)
        self.play(FadeOut(VGroup(heading, input_card, input_label, cnn_group, moving_heatmap, prob_title, probabilities, comparison, match_group, architecture_group)), run_time=0.8)

    def closing_scene(self):
        heading = self._label("Pipeline overview and test performance", 34, WHITE, BOLD).to_edge(UP, buff=0.38)
        steps = [
            ("1000G VCFs", BLUE),
            ("H1 windows", CYAN),
            ("AlphaGenome", PURPLE),
            ("32 kbp heatmaps", GREEN),
            ("per-track norm", YELLOW),
            ("CNN", PINK),
            ("Superpopulation", ORANGE),
        ]
        nodes = VGroup()
        for i, (text, color) in enumerate(steps):
            dot = Circle(radius=0.055, stroke_width=0, fill_color=color, fill_opacity=0.90)
            label = Text(text, font_size=14, color=WHITE if i in [0, 2, 5, 6] else MUTED, weight=BOLD if i in [0, 2, 5, 6] else NORMAL)
            node = VGroup(dot, label).arrange(DOWN, buff=0.08)
            nodes.add(node)
        nodes.arrange(RIGHT, buff=0.42).move_to([0, 1.95, 0])
        rail = Line(nodes[0][0].get_center(), nodes[-1][0].get_center(), color=MUTED, stroke_width=2.0, stroke_opacity=0.35)
        pipeline_group = VGroup(rail, nodes)

        card = RoundedRectangle(corner_radius=0.18, width=8.00, height=3.30, stroke_color=CYAN, stroke_width=1.4, fill_color=CYAN, fill_opacity=0.045).move_to([0, -0.55, 0])
        result_title = Text("Classification performance on the test set", font_size=25, color=WHITE, weight=BOLD).move_to([0, 0.74, 0])
        result_note = Text('"W-" denotes weighted-average.', font_size=16, color=MUTED).next_to(result_title, DOWN, buff=0.12)
        metrics = [
            ("W-Precision", "0.75", CYAN),
            ("W-Recall", "0.74", GREEN),
            ("W-F1-score", "0.74", YELLOW),
            ("Accuracy", "0.74", PINK),
        ]
        table = VGroup()
        for metric, value, color in metrics:
            value_text = Text(value, font_size=34, color=color, weight=BOLD)
            metric_text = Text(metric, font_size=16, color=WHITE, weight=BOLD).next_to(value_text, DOWN, buff=0.10)
            cell = VGroup(value_text, metric_text)
            halo = RoundedRectangle(corner_radius=0.12, width=1.28, height=0.96, stroke_color=color, stroke_width=0.8, fill_color=color, fill_opacity=0.075).move_to(cell)
            table.add(VGroup(halo, cell))
        table.arrange(RIGHT, buff=0.22).move_to([0, -0.88, 0])
        result_group = VGroup(card, result_title, result_note, table)
        self.play(Write(heading), Create(rail), LaggedStart(*[FadeIn(n, scale=0.92) for n in nodes], lag_ratio=0.08), run_time=1.6)
        self.play(FadeIn(result_group, shift=UP * 0.18), run_time=1.4)
        self.wait(11.0)
