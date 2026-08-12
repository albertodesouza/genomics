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
source scripts/env/start_genomics_universal.sh
scripts/env/install_genomics_env.sh
source scripts/maintenance/vep_install.sh
scripts/ops/run_in_background.sh --config configs/genomes_analyzer/config_human_30x_latest_ref.yaml
scripts/ops/monitor_monster.sh
scripts/diagnostics/diagnose_bcftools_error.sh
```
