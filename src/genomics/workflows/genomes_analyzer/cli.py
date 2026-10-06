"""Command-line interface for the genomes analyzer pipeline."""

import argparse

import yaml
from rich.panel import Panel

from . import legacy
from . import pipeline
from .state import set_config


def cli_main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Pipeline humano — Low-Memory: índice BWA pré-pronto, verbose, "
            "cache, CRAM, cobertura e downsample"
        )
    )
    parser.add_argument("-c", "--config", required=True, help="Arquivo YAML de configuração")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        cfg_raw = yaml.safe_load(f)

    cfg = pipeline.normalize_config_schema(cfg_raw)
    set_config(cfg)
    legacy.console.print(
        Panel.fit(
            "[bold]YAML normalizado[/bold]\n"
            f"Chaves topo: {list(cfg.keys())}\n"
            f"General → {', '.join(sorted(list(cfg['general'].keys()))[:12])}...",
            border_style="cyan",
        )
    )
    pipeline.main(cfg)


if __name__ == "__main__":
    cli_main()
