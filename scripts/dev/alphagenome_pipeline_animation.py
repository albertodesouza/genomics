from pathlib import Path

import numpy as np
from manim import *


ASSET_DIR = Path(__file__).resolve().parents[2] / "assets" / "pipeline_animation" / "ppt_media"

BG = "#070B1A"
CYAN = "#41D9FF"
BLUE = "#4E7BFF"
PURPLE = "#9B5CFF"
PINK = "#FF4FA3"
GREEN = "#59E39B"
YELLOW = "#FFD166"
ORANGE = "#FF9F1C"
MUTED = "#8EA0C5"


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

    def _dna(self, width=5.2, bases=18):
        pair_groups = VGroup()
        xs = np.linspace(-width / 2, width / 2, bases)
        pairs = [("A", "T"), ("C", "G"), ("G", "C"), ("T", "A")]
        top_points = []
        bot_points = []
        for i, x in enumerate(xs):
            y = 0.28 * np.sin(i * 0.9)
            top_points.append([x, y, 0])
            bot_points.append([x, -y, 0])
            rung = Line([x, y, 0], [x, -y, 0], color=BLUE, stroke_width=1.5, stroke_opacity=0.65)
            top_base, bottom_base = pairs[i % len(pairs)]
            top_offset = 0.16 if y >= 0 else -0.16
            bottom_offset = -0.16 if y >= 0 else 0.16
            top_label = Text(top_base, font_size=13, color=CYAN if i % 5 else PINK).move_to([x, y + top_offset, 0])
            bottom_label = Text(bottom_base, font_size=13, color=PURPLE if i % 5 else PINK).move_to([x, -y + bottom_offset, 0])
            pair_groups.add(VGroup(rung, top_label, bottom_label))
        top = VMobject(color=CYAN, stroke_width=3)
        bottom = VMobject(color=PURPLE, stroke_width=3)
        top.set_points_smoothly(top_points)
        bottom.set_points_smoothly(bot_points)
        group = VGroup(top, bottom, pair_groups)
        group.strands = VGroup(top, bottom)
        group.base_pairs = pair_groups
        return group

    def _dna_growth(self, dna, run_time=1.5):
        return AnimationGroup(
            Create(dna.strands, lag_ratio=0.0),
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
        strands = VGroup()
        strand_colors = [CYAN, GREEN, YELLOW, PINK, PURPLE, ORANGE, BLUE, CYAN, GREEN]
        for i, color in enumerate(strand_colors):
            strand = self._single_strand(width=4.15, bases=16, color=color).scale(0.42)
            strand.shift(RIGHT * 1.35 + UP * (1.35 - i * 0.32))
            strands.add(strand)
        h1 = self._badge("H1 haplotype sequences", PINK, width=2.8).next_to(strands, DOWN, buff=0.34)
        brace = Brace(strands, RIGHT, color=YELLOW)
        count = Text("3,202 individuals", font_size=24, color=YELLOW, weight=BOLD).next_to(brace, RIGHT, buff=0.16)
        arrows = VGroup(self._flow_arrow(vcf.get_right(), strands.get_left() + UP * 0.45), self._flow_arrow(ref.get_right(), strands.get_left() + DOWN * 0.45, MUTED))
        self.play(Write(heading), FadeIn(vcf), FadeIn(ref), run_time=1.0)
        self.play(Create(arrows), run_time=0.9)
        self.play(LaggedStart(*[self._single_strand_growth(strand, run_time=0.75) for strand in strands], lag_ratio=0.11), FadeIn(h1), run_time=1.7)
        self.play(Create(brace), FadeIn(count), run_time=0.8)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, vcf, ref, arrows, strands, h1, brace, count)), run_time=0.7)

    def alphagenome_scene(self):
        heading = self._label("AlphaGenome: one strand in, 5930 human tracks out", 32, WHITE, BOLD).to_edge(UP, buff=0.45)
        helix = self._dna(width=2.9, bases=12).scale(0.58).shift(LEFT * 4.65 + UP * 1.05)
        plus = self._single_strand(width=3.1, bases=13, color=CYAN).scale(0.82).move_to(helix).shift(DOWN * 1.25)
        plus_label = self._badge("+ strand pass", BLUE, width=1.65).next_to(plus, DOWN, buff=0.12)
        minus_label = self._badge("- strand later", PURPLE, width=1.6).next_to(plus_label, DOWN, buff=0.16)
        model = RoundedRectangle(corner_radius=0.18, width=2.55, height=1.35, stroke_color=CYAN, fill_color=PURPLE, fill_opacity=0.18)
        model_text = VGroup(Text("AlphaGenome", font_size=27, color=WHITE, weight=BOLD), Text("single-strand pass", font_size=16, color=MUTED)).arrange(DOWN, buff=0.12)
        model_text.move_to(model)
        model_group = VGroup(model, model_text)
        model_group.shift(RIGHT * 0.25)
        tracks = VGroup()
        track_colors = [CYAN, GREEN, YELLOW, PINK, PURPLE, ORANGE, BLUE]
        for i, color in enumerate(track_colors):
            tr = self._signal_track(width=2.95, height=0.34, color=color, phase=i).scale(0.74).shift(RIGHT * 4.55 + UP * (1.7 - i * 0.38))
            tracks.add(tr)
        count = self._badge("5930 human tracks", ORANGE, width=2.35).next_to(tracks, DOWN, buff=0.22)
        ellipsis = Text("...", font_size=34, color=MUTED).next_to(tracks[-1], DOWN, buff=0.02)
        arrow_in = self._flow_arrow(plus.get_right(), model_group.get_left(), BLUE)
        arrow_out = self._flow_arrow(model_group.get_right(), tracks.get_left(), CYAN)
        self.play(Write(heading), self._dna_growth(helix), run_time=1.0)
        self.play(ReplacementTransform(helix.strands[0].copy(), plus.strand), LaggedStart(*[FadeIn(b) for b in plus.bases], lag_ratio=0.04), FadeIn(plus_label), FadeIn(minus_label), run_time=1.3)
        aligned = plus.copy().scale(0.78).next_to(model_group, LEFT, buff=0.35)
        self.play(FadeIn(model_group, scale=0.96), Transform(plus, aligned), plus_label.animate.next_to(aligned, DOWN, buff=0.12), run_time=1.0)
        self.play(Create(arrow_in), plus.animate.move_to(model_group), plus.animate.set_opacity(0.0), run_time=1.0)
        self.play(Create(arrow_out), LaggedStart(*[Create(t) for t in tracks], lag_ratio=0.08), FadeIn(ellipsis), FadeIn(count, shift=UP * 0.15), run_time=1.8)
        self.wait(1.8)
        self.play(FadeOut(VGroup(heading, helix, plus, plus_label, minus_label, model_group, tracks, ellipsis, count, arrow_in, arrow_out)), run_time=0.7)

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
