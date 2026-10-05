# Environment

The `genomics` environment is created by `scripts/env/install.sh` (see [Installation](installation.md)).
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

## AlphaGenome server environment

The local AlphaGenome model server does not run in `genomics`: it needs JAX with CUDA, TensorFlow and
`alphagenome>=0.7`, which conflict with PyTorch. `genomics alphagenome server setup` creates a second
conda environment, `alphagenome` (or a venv inside the `alphagenome_research` checkout when there is no
conda), and the visualizer and `genomics alphagenome server start` run the server with its interpreter.
You never need to activate it. Override the choice with `ALPHAGENOME_SERVER_PYTHON` /
`ALPHAGENOME_SERVER_DIR`.

| Variable | Default | Meaning |
|---|---|---|
| `ALPHAGENOME_API_KEY` | – (also read from `~/.env`) | hosted API key |
| `ALPHAGENOME_ADDRESS` | – | use a self-hosted server, e.g. `grpc://127.0.0.1:50051` |
| `ALPHAGENOME_TLS_CA_CERT` | – | CA certificate of a self-signed TLS server |
| `ALPHAGENOME_SERVER_DIR` | `../alphagenome_research`, then `~/.local/share/genomics/alphagenome_research` | checkout used for the local server |
| `ALPHAGENOME_SERVER_PYTHON` | the env `alphagenome` or another env with CUDA `jax` | interpreter of the local server |
| `ALPHAGENOME_SERVER_HOST` / `ALPHAGENOME_SERVER_PORT` | `127.0.0.1` / `50051` | where the local server listens |
| `GENOMICS_DATA_ROOT` / `GENOMICS_RESULTS_ROOT` | `/dados/GENOMICS_DATA` / `results/` | datasets / runs and caches |
