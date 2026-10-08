# Scripts

Root-level script wrappers have been removed. Use categorized paths directly.

| Category | Path | Purpose |
|---|---|---|
| Environment | `scripts/env/` | Conda/bootstrap/install scripts |
| Operations | `scripts/ops/` | background runs, monitors, production helpers |
| Maintenance | `scripts/maintenance/` | VEP/reference/dependency maintenance |
| Diagnostics | `scripts/diagnostics/` | structure/API checks and debug probes |
| Experiments | `scripts/experiments/` | benchmark and overnight experiment runners |
| Development | `scripts/dev/` | demos, plotting checks, local tests |

## Common Scripts

```bash
scripts/env/install.sh                      # visualizer env (bcftools, samtools, [visualizer]); --training, --alphagenome-server, --no-conda
source scripts/env/start_genomics_universal.sh
scripts/env/install_genomics_env.sh         # full bioinformatics toolchain for genomes-analyzer
source scripts/maintenance/vep_install.sh
scripts/ops/run_in_background.sh --config configs/genomes_analyzer/config_human_30x_latest_ref.yaml
scripts/ops/monitor_monster.sh
scripts/diagnostics/diagnose_bcftools_error.sh
python3 scripts/dev/capture_docs_screenshots.py --url http://127.0.0.1:8780/   # retake the docs screenshots from a running visualizer
```

`capture_docs_screenshots.py` drives the visualizer in headless Chromium (Playwright) and writes the
WebP screenshots and figures under `docs/assets/visualizer/`. The shots follow the tutorial's
running example on `1kg_high_coverage`. Some need an AlphaGenome backend (Perturbation Lab,
expression in other tissues), a loadable trained run, and network access for GTEx, Ensembl and
ENCODE. `--list` shows the shots and `--only a,b` retakes some.
