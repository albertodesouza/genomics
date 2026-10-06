"""Allow running the pipeline with `python -m genomics.workflows.genomes_analyzer`."""

from .cli import cli_main


if __name__ == "__main__":
    cli_main()
