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
        self.vcf_scene()
        self.gene_window_scene()
        self.reconstruction_scene()
        self.alphagenome_scene()
        self.track_selection_scene()
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
            labels.add(Text(letter, font_size=font_size, color=label_color).move_to(np.array(point) + letter_offset))
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
        top = self._strand_from_points(top_points, top_letters, color=CYAN, letter_offset=UP * 0.16)
        bottom = self._strand_from_points(bot_points, bottom_letters, color=PURPLE, letter_offset=DOWN * 0.16)
        group = VGroup(top, bottom, pair_groups)
        group.top_strand = top
        group.bottom_strand = bottom
        group.strands = VGroup(top, bottom)
        group.base_pairs = pair_groups
        return group

    def _dna_growth(self, dna, run_time=1.5):
        return AnimationGroup(
            Create(VGroup(dna.top_strand.curve, dna.bottom_strand.curve), lag_ratio=0.0),
            LaggedStart(*[FadeIn(base, scale=0.9) for base in VGroup(dna.top_strand.bases, dna.bottom_strand.bases)], lag_ratio=0.035),
            LaggedStart(*[FadeIn(pair, scale=0.9) for pair in dna.base_pairs], lag_ratio=0.045),
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

    def title_scene(self):
        title = Text("AlphaGenome-Derived Functional Signals", font_size=42, color=WHITE, weight=BOLD)
        subtitle = Text("Inferring superpopulation labels from individual DNA variation", font_size=24, color=MUTED)
        dna = self._dna(width=7.0, bases=28).scale(0.9).shift(DOWN * 1.0)
        group = VGroup(title, subtitle).arrange(DOWN, buff=0.22).shift(UP * 1.35)
        self.play(FadeIn(group, shift=UP * 0.25), self._dna_growth(dna), run_time=2.2)
        self.wait(2.5)
        self.play(FadeOut(group), FadeOut(dna), run_time=0.8)

    def vcf_scene(self):
        heading = self._label("1000 Genomes high-coverage input", 34, WHITE, BOLD).to_edge(UP, buff=0.45)
        cards = VGroup()
        for i, label in enumerate(["chr1.vcf.gz", "chr2.vcf.gz", "chr3.vcf.gz", "...", "chr22.vcf.gz"]):
            rect = RoundedRectangle(corner_radius=0.12, width=1.65, height=0.72, stroke_color=CYAN, fill_color=BLUE, fill_opacity=0.12)
            txt = Text(label, font_size=18, color=WHITE)
            txt.move_to(rect)
            card = VGroup(rect, txt)
            card.move_to([i * 1.75 - 3.5, 1.25, 0])
            cards.add(card)
        vcf_example = self._maybe_image("image18.png", height=2.45, opacity=0.72).shift(DOWN * 1.0)
        example_frame = RoundedRectangle(
            corner_radius=0.10,
            width=5.9,
            height=2.65,
            stroke_color=CYAN,
            stroke_width=1.2,
            stroke_opacity=0.70,
            fill_color=BLUE,
            fill_opacity=0.06,
        ).move_to(vcf_example)
        example_panel = Group(example_frame, vcf_example)
        self.play(Write(heading), LaggedStart(*[FadeIn(c, shift=DOWN * 0.2) for c in cards], lag_ratio=0.12), run_time=1.6)
        self.play(FadeIn(example_frame), FadeIn(vcf_example), run_time=1.1)
        self.wait(2.0)
        self.play(FadeOut(Group(heading, cards, example_panel)), run_time=0.6)

    def gene_window_scene(self):
        heading = self._label("Center 512 kbp windows on pigmentation genes", 33, WHITE, BOLD).to_edge(UP, buff=0.45)
        chromosome = Line(LEFT * 5.8, RIGHT * 5.8, color=MUTED, stroke_width=7, stroke_opacity=0.55).shift(UP * 0.25)
        genes = ["DDB1", "EDAR", "HERC2", "MC1R", "MFSD12", "OCA2", "SLC24A5", "SLC45A2", "TCHH", "TYR", "TYRP1"]
        markers = VGroup()
        for i, gene in enumerate(genes):
            x = -5.25 + i * 1.05
            color = [CYAN, BLUE, PURPLE, PINK, GREEN, YELLOW][i % 6]
            tick = Line([x, 0.0, 0], [x, 0.65, 0], color=color, stroke_width=4)
            txt = Text(gene, font_size=17, color=color).rotate(PI / 5).next_to(tick, UP, buff=0.06)
            markers.add(VGroup(tick, txt))
        bracket = VGroup(
            Line([-0.38, -0.55, 0], [0.38, -0.55, 0], color=YELLOW, stroke_width=4),
            Line([-0.38, -0.55, 0], [-0.38, -0.18, 0], color=YELLOW, stroke_width=4),
            Line([0.38, -0.55, 0], [0.38, -0.18, 0], color=YELLOW, stroke_width=4),
            Text("512 kbp", font_size=23, color=YELLOW).shift(DOWN * 0.95),
        )
        self.play(Write(heading), Create(chromosome), run_time=1.1)
        self.play(LaggedStart(*[FadeIn(m, shift=UP * 0.15) for m in markers], lag_ratio=0.06), run_time=1.7)
        self.play(Create(bracket), run_time=0.8)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, chromosome, markers, bracket)), run_time=0.7)

    def reconstruction_scene(self):
        heading = self._label("Reconstruct individual H1 haplotype sequences", 33, WHITE, BOLD).to_edge(UP, buff=0.45)
        vcf = self._badge("VCF variants", BLUE, width=2.0).shift(LEFT * 4.25 + UP * 0.55)
        ref = self._badge("reference genome", MUTED, width=2.35).shift(LEFT * 4.25 + DOWN * 0.45)
        merge_point = LEFT * 1.75
        h1_sequences = VGroup()
        helix_colors = [CYAN, GREEN, YELLOW, PINK, PURPLE, ORANGE, BLUE, CYAN, GREEN]
        for i, color in enumerate(helix_colors):
            helix = self._dna(width=2.25, bases=10).scale(0.34)
            helix.top_strand.curve.set_color(color)
            helix.bottom_strand.curve.set_color(color)
            row, col = divmod(i, 3)
            helix.move_to([0.25 + col * 1.85, 1.15 - row * 0.78, 0])
            h1_sequences.add(helix)
        h1 = self._badge("H1 haplotype sequences", PINK, width=2.8).next_to(h1_sequences, DOWN, buff=0.34)
        count = Text("3,202 individual H1 sequences", font_size=24, color=YELLOW, weight=BOLD).next_to(h1, DOWN, buff=0.20)
        self.play(Write(heading), FadeIn(vcf), FadeIn(ref), run_time=1.0)
        self.play(vcf.animate.move_to(merge_point + UP * 0.20), ref.animate.move_to(merge_point + DOWN * 0.20), run_time=0.9)
        self.play(vcf.animate.scale(0.78).set_opacity(0.0), ref.animate.scale(0.78).set_opacity(0.0), run_time=0.45)
        self.play(LaggedStart(*[self._dna_growth(helix, run_time=0.8) for helix in h1_sequences], lag_ratio=0.09), FadeIn(h1), run_time=1.9)
        self.play(FadeIn(count, shift=UP * 0.08), run_time=0.55)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, vcf, ref, h1_sequences, h1, count)), run_time=0.7)

    def alphagenome_scene(self):
        heading = self._label("AlphaGenome: one strand in, 5930 human tracks out", 32, WHITE, BOLD).to_edge(UP, buff=0.45)
        helix = self._dna(width=2.9, bases=12).scale(0.58).shift(LEFT * 4.95)
        plus_sequence = "ACGTACGTACGT"
        minus_sequence = "TGCATGCATGCA"
        plus = self._linear_strand(plus_sequence, UP * 0.95, DOWN * 0.95, color=CYAN).shift(LEFT * 1.85 + UP * 1.25)
        minus = self._linear_strand(minus_sequence, UP * 0.95, DOWN * 0.95, color=PURPLE).shift(LEFT * 1.85 + DOWN * 1.25)
        models = VGroup()
        for y in [1.25, -1.25]:
            model = RoundedRectangle(corner_radius=0.18, width=2.55, height=1.1, stroke_color=CYAN, fill_color=PURPLE, fill_opacity=0.18)
            model_text = VGroup(Text("AlphaGenome", font_size=25, color=WHITE, weight=BOLD), Text("strand model", font_size=15, color=MUTED)).arrange(DOWN, buff=0.10)
            model_text.move_to(model)
            model_group = VGroup(model, model_text).shift(RIGHT * 0.45 + UP * y)
            models.add(model_group)
        top_tracks = VGroup()
        bottom_tracks = VGroup()
        track_colors = [CYAN, GREEN, YELLOW, PINK, PURPLE, ORANGE, BLUE]
        for i, color in enumerate(track_colors):
            top_tracks.add(self._signal_track(width=2.35, height=0.24, color=color, phase=i).scale(0.58).shift(RIGHT * 4.65 + UP * (2.05 - i * 0.19)))
            bottom_tracks.add(self._signal_track(width=2.35, height=0.24, color=color, phase=i + 0.7).scale(0.58).shift(RIGHT * 4.65 + DOWN * (0.45 + i * 0.19)))
        top_count = self._badge("5930 tracks", ORANGE, width=1.62).next_to(top_tracks, DOWN, buff=0.14)
        bottom_count = self._badge("5930 tracks", ORANGE, width=1.62).next_to(bottom_tracks, DOWN, buff=0.14)
        ellipses = VGroup(Text("...", font_size=26, color=MUTED).next_to(top_tracks[-1], DOWN, buff=0.01), Text("...", font_size=26, color=MUTED).next_to(bottom_tracks[-1], DOWN, buff=0.01))
        self.play(Write(heading), self._dna_growth(helix), run_time=1.0)
        self.play(
            FadeOut(helix.base_pairs),
            helix.top_strand.animate.shift(UP * 0.24),
            helix.bottom_strand.animate.shift(DOWN * 0.24),
            run_time=0.55,
        )
        self.play(
            Transform(helix.top_strand, plus),
            Transform(helix.bottom_strand, minus),
            run_time=1.3,
        )
        plus = helix.top_strand
        minus = helix.bottom_strand
        self.play(FadeIn(models, scale=0.96), run_time=0.75)
        self.play(plus.animate.move_to(models[0]).set_opacity(0.0), minus.animate.move_to(models[1]).set_opacity(0.0), run_time=1.0)
        self.play(
            LaggedStart(*[Create(t) for t in top_tracks], lag_ratio=0.08),
            LaggedStart(*[Create(t) for t in bottom_tracks], lag_ratio=0.08),
            FadeIn(ellipses),
            FadeIn(top_count, shift=UP * 0.12),
            FadeIn(bottom_count, shift=UP * 0.12),
            run_time=1.8,
        )
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, helix, plus, minus, models, top_tracks, bottom_tracks, ellipses, top_count, bottom_count)), run_time=0.7)

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
        strands = VGroup(self._badge("+ strand", BLUE, width=1.4), self._badge("- strand", PURPLE, width=1.4)).arrange(RIGHT, buff=0.45).next_to(labels, DOWN, buff=0.3)
        self.play(Write(heading), LaggedStart(*[Create(t) for t in tracks], lag_ratio=0.06), run_time=1.5)
        self.play(FadeIn(source, shift=UP * 0.2), run_time=0.6)
        self.play(tracks.animate.scale(0.92), source.animate.set_opacity(0.22), FadeIn(labels, shift=UP * 0.2), FadeIn(strands, shift=UP * 0.2), run_time=1.3)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, tracks, source, labels, strands)), run_time=0.7)

    def normalization_scene(self):
        heading = self._label("Crop the center 32 kbp and log-normalize", 33, WHITE, BOLD).to_edge(UP, buff=0.45)
        raw = self._signal_track(width=8.5, color=CYAN, phase=0.4).shift(UP * 1.0)
        crop = Rectangle(width=2.0, height=1.05, stroke_color=YELLOW, stroke_width=4, fill_color=YELLOW, fill_opacity=0.06).move_to(raw)
        crop_label = Text("centered 32 kbp", font_size=23, color=YELLOW).next_to(crop, DOWN, buff=0.22)
        formula = MathTex(r"x_{norm}=\frac{\ln(x+1)}{\ln(x_{max}+1)}", font_size=42, color=WHITE).shift(DOWN * 1.0)
        normalized = self._signal_track(width=4.2, color=GREEN, phase=0.9).shift(DOWN * 2.15)
        norm_label = Text("normalized to comparable signal scale", font_size=22, color=GREEN).next_to(normalized, DOWN, buff=0.15)
        self.play(Write(heading), Create(raw), run_time=1.0)
        self.play(Create(crop), FadeIn(crop_label), run_time=0.9)
        self.play(Write(formula), run_time=1.1)
        self.play(ReplacementTransform(crop.copy(), normalized), FadeIn(norm_label), run_time=1.2)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, raw, crop, crop_label, formula, normalized, norm_label)), run_time=0.7)

    def dataset_scene(self):
        heading = self._label("Concatenate signals into a supervised dataset", 33, WHITE, BOLD).to_edge(UP, buff=0.45)
        heatmaps = VGroup()
        labels = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        colors = [PINK, ORANGE, YELLOW, CYAN, GREEN]
        for i in range(5):
            hm = self._heatmap(rows=8, cols=14, cell=0.115)
            hm.shift(LEFT * 4.2 + RIGHT * i * 2.1 + UP * (0.4 if i % 2 else 0.05))
            badge = self._badge(labels[i], colors[i], width=0.86).next_to(hm, DOWN, buff=0.15)
            heatmaps.add(VGroup(hm, badge))
        stack_label = Text("genes x tissues x strands x positions", font_size=24, color=MUTED).shift(DOWN * 2.15)
        self.play(Write(heading), LaggedStart(*[FadeIn(h, shift=UP * 0.25) for h in heatmaps], lag_ratio=0.12), run_time=1.8)
        self.play(FadeIn(stack_label), heatmaps.animate.arrange(RIGHT, buff=0.2).scale(0.9).shift(UP * 0.15), run_time=1.2)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, heatmaps, stack_label)), run_time=0.7)

    def closing_scene(self):
        heading = self._label("Functional-signal pipeline", 34, WHITE, BOLD).to_edge(UP, buff=0.4)
        steps = [
            ("1000G VCFs", BLUE),
            ("H1 windows", CYAN),
            ("AlphaGenome", PURPLE),
            ("RNA-seq tracks", GREEN),
            ("32 kbp norm", YELLOW),
            ("Classifier", PINK),
            ("Superpopulation", ORANGE),
        ]
        nodes = VGroup(*[self._badge(text, color, width=1.55 if len(text) < 12 else 1.9) for text, color in steps]).arrange(RIGHT, buff=0.2).scale(0.68)
        nodes.shift(UP * 0.65)
        arrows = VGroup()
        for left, right in zip(nodes[:-1], nodes[1:]):
            arrows.add(self._flow_arrow(left.get_right(), right.get_left(), CYAN).scale(0.8))
        classifier = VGroup()
        for i, label in enumerate(["AFR", "AMR", "EAS", "EUR", "SAS"]):
            bar = Rectangle(width=0.35, height=0.35 + 0.22 * ((i * 3) % 5), fill_color=[PINK, ORANGE, YELLOW, CYAN, GREEN][i], fill_opacity=0.8, stroke_width=0)
            txt = Text(label, font_size=16, color=WHITE).next_to(bar, DOWN, buff=0.08)
            classifier.add(VGroup(bar, txt))
        classifier.arrange(RIGHT, buff=0.22).shift(DOWN * 1.45)
        message = Text("DNA variation transformed into predicted function, then learned as ancestry-informative signal", font_size=23, color=WHITE).to_edge(DOWN, buff=0.45)
        self.play(Write(heading), LaggedStart(*[FadeIn(n, scale=0.92) for n in nodes], lag_ratio=0.08), run_time=1.5)
        self.play(LaggedStart(*[Create(a) for a in arrows], lag_ratio=0.08), run_time=0.9)
        self.play(FadeIn(classifier, shift=UP * 0.25), Write(message), run_time=1.3)
        self.wait(8.5)
