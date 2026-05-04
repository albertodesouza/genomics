# Views

`views/` guarda definicoes reutilizaveis de subconjuntos logicos do dataset fisico canonico.

Uma view deve descrever apenas a selecao e a forma de entrada dos dados, por exemplo:

- `dataset_dir`
- `sample_ids`
- `sample_ids_path`
- `superpopulations_to_use`
- `populations_to_use`
- `genes_to_use`
- `alphagenome_outputs`
- `haplotype_mode`
- `window_center_size`
- `downsample_factor`
- `normalization_method`
- `tensor_layout`
- `ontology_terms`
- `indel_include_valid_mask`
- `indel_neutral_value`

Uma view nao deve carregar decisoes de experimento, por exemplo:

- `prediction_target`
- `known_classes`
- `derived_targets`
- `model`
- `training`
- `data_split`
- `wandb`

Esses itens continuam em `configs/`.

## Regra pratica

- `views/`: o que entra no modelo
- `configs/`: o que o modelo faz com isso

## Uso

No YAML:

```yaml
dataset_input:
  view_path: "/abs/path/to/my.view.json"
  processed_cache_dir: "/dados/.../runs"
```

Os campos definidos no YAML sobrescrevem os defaults da view.
