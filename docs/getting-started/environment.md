# Environment

Use the project activation script instead of calling `conda activate genomics` manually:

```bash
source scripts/env/start_genomics_universal.sh
```

The universal script:

- detects common Conda/Miniforge locations;
- activates the `genomics` environment;
- cleans problematic CUDA library paths on affected systems;
- loads Bash completion for `genomics` when the command is installed.

If your Conda installation is always under `~/miniforge3`, this simpler script is also available:

```bash
source scripts/env/start_genomics.sh
```

## Completion

The activation scripts run this automatically:

```bash
source <(genomics completion bash)
```

To install completion persistently for your user:

```bash
mkdir -p ~/.local/share/bash-completion/completions
genomics completion bash > ~/.local/share/bash-completion/completions/genomics
```
