# Scripts Layout

Scripts are organized by category. Use the categorized paths directly; root-level compatibility wrappers have been removed.

| Category | Intended path | Examples |
|---|---|---|
| Environment setup | `scripts/env/` | Conda/bootstrap/install scripts |
| Operations | `scripts/ops/` | background runs, monitors, production helpers |
| Maintenance | `scripts/maintenance/` | reference fixes, VEP updates, data repair |
| Diagnostics | `scripts/diagnostics/` | structure/API checks and debug probes |
| Experiments | `scripts/experiments/` | benchmark and overnight experiment runners |
| Development | `scripts/dev/` | plotting tests, demos, local checks |

## Current Mapping

| Category | Scripts |
|---|---|
| `env/` | `start_genomics_universal.sh`, `start_genomics.sh`, `install_genomics_env.sh`, `install_genomics_env_linux_aarch64_nvidia.sh`, `install_conda_universal.sh`, `conda_initialize.sh`, `install_alphagenome.sh` |
| `ops/` | `run_atena_background.sh`, `run_monster_background.sh`, `run_in_background.sh`, `monitor_monster.sh`, `monitor_bwa_index.sh`, `setup_monster_complete.sh`, `setup_monster_256gb.sh`, `pipeline.sh`, `run.bat` |
| `maintenance/` | `vep_install_fixed.sh`, `vep_install_smart.sh`, `vep_install.sh`, `vep_install_latest.sh`, `vep_update_to_latest.sh`, `fix_reference_mismatch.sh`, `fix_vep_dependencies.sh` |
| `diagnostics/` | `check_tsv_structure.py`, `diagnose_bcftools_error.sh`, `check_output_structure.py`, `pairwise_vcf_diff.sh`, `read_alphagenome_predictions.py`, `check_alphagenome_outputs.py`, `check_visualization_api.py`, `check_dna_client_methods.py` |
| `experiments/` | `run_important_experiments_overnight.sh`, `test_parallel.sh` |
| `dev/` | `demo_csv_feature.sh`, `test_simple_plot.py`, `test_manual_plot.py` |

When updating automation, use the categorized paths shown above.
