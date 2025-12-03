# Labels de Ontologias Dinâmicos

## Mudança Implementada

Modificado `verify_processed_dataset.py` para gerar labels de ontologias dinamicamente a partir dos metadados do dataset, eliminando a necessidade de manter labels hardcoded.

## Problema Anterior

**Antes**: Labels de ontologias estavam hardcoded:

```python
ONTOLOGY_LABELS = [
    "CL:1000458 (+)\nMelanocyte",
    "CL:0000346 (+)\nDermal Papilla",
    "CL:2000092 (+)\nKeratinocyte",
    "CL:1000458 (-)\nMelanocyte",
    "CL:0000346 (-)\nDermal Papilla",
    "CL:2000092 (-)\nKeratinocyte"
]
```

**Limitações**:
- Labels precisavam ser atualizados manualmente para cada dataset
- Nomes simplificados não refletiam os nomes completos dos biosamples
- Inconsistente se dataset tivesse ontologias diferentes

## Solução Implementada

**Agora**: Geração dinâmica de labels a partir dos metadados do dataset:

```python
def generate_ontology_labels(dataset_metadata: Dict) -> List[str]:
    """
    Gera labels de ontologias dinamicamente a partir dos metadados do dataset.
    
    Lê o campo 'ontology_details' que contém informações completas sobre cada
    ontologia, incluindo biosample_name, biosample_type, etc.
    """
    ontology_details = dataset_metadata.get('ontology_details', {})
    
    if not ontology_details:
        # Fallback para labels padrão se metadados não tiverem ontology_details
        return [labels hardcoded como fallback]
    
    sorted_ontologies = sorted(ontology_details.keys())
    labels = []
    
    # Primeiro todas as ontologias com strand +
    for ontology_curie in sorted_ontologies:
        details = ontology_details[ontology_curie]
        biosample_name = details.get('biosample_name', ontology_curie)
        label = f"{ontology_curie} (+)\n{biosample_name}"
        labels.append(label)
    
    # Depois todas as ontologias com strand -
    for ontology_curie in sorted_ontologies:
        details = ontology_details[ontology_curie]
        biosample_name = details.get('biosample_name', ontology_curie)
        label = f"{ontology_curie} (-)\n{biosample_name}"
        labels.append(label)
    
    return labels
```

## Fonte de Dados

Os labels são gerados a partir do campo `ontology_details` nos metadados do dataset:

```json
"ontology_details": {
    "CL:0000346": {
        "biosample_name": "hair follicle dermal papilla cell",
        "biosample_type": "primary_cell",
        "biosample_life_stage": "adult",
        ...
    },
    "CL:1000458": {
        "biosample_name": "melanocyte of skin",
        "biosample_type": "primary_cell",
        ...
    },
    ...
}
```

## Resultado

Labels gerados automaticamente:

```
1. CL:0000346 (+) | hair follicle dermal papilla cell
2. CL:1000458 (+) | melanocyte of skin
3. CL:2000092 (+) | hair follicular keratinocyte
4. CL:0000346 (-) | hair follicle dermal papilla cell
5. CL:1000458 (-) | melanocyte of skin
6. CL:2000092 (-) | hair follicular keratinocyte
```

## Benefícios

1. **Dinâmico**: Labels refletem automaticamente as ontologias do dataset
2. **Nomes completos**: Usa biosample_name oficial em vez de apelidos
3. **Compatibilidade**: Mantém fallback para datasets antigos sem ontology_details
4. **Manutenibilidade**: Não requer edição de código ao mudar ontologias
5. **Consistência**: Labels sempre sincronizados com o dataset real
6. **Flexibilidade**: Funciona com qualquer conjunto de ontologias

## Funções Modificadas

Três funções de plotagem foram modificadas para aceitar `dataset_metadata`:

1. **`plot_comparison()`** - Adiciona comparação de tracks cache vs AlphaGenome
2. **`plot_individual_comparison()`** - Compara dois indivíduos
3. **`plot_raw_data()`** - Mostra dados brutos do AlphaGenome

Todas agora têm parâmetro opcional `dataset_metadata: Optional[Dict] = None`.

## Integração

A função `generate_ontology_labels()` é chamada dentro de cada função de plot:

```python
# Gerar labels dinamicamente dos metadados
if dataset_metadata:
    ontology_labels = generate_ontology_labels(dataset_metadata)
else:
    ontology_labels = generate_ontology_labels({})  # Usa fallback
```

O `dataset_metadata` é obtido automaticamente do viewer quando disponível:

```python
dataset_metadata = getattr(viewer, 'dataset_metadata', None) if viewer else None
```

## Compatibilidade

✅ **Retrocompatível**: Fallback para labels hardcoded se campo não existir
✅ **Sem quebras**: Funciona com datasets antigos e novos
✅ **Validação**: Exibe warning se usar fallback

## Teste de Verificação

```bash
✅ Ontology details encontradas nos metadados: 3 ontologias
✅ Labels gerados dinamicamente: 6 labels
✅ Função generate_ontology_labels funcionará corretamente!
```

## Arquivo Modificado

- `neural_ancestry_predictor/verify_processed_dataset.py`

## Requer Metadados Atualizados

Para aproveitar essa funcionalidade, o dataset deve ter sido gerado com a versão atualizada do `dataset_builder.py` que inclui o campo `ontology_details` nos metadados.

Para datasets antigos sem este campo, a função usa fallback com labels padrão e exibe um warning.

---

## Update: Redução do Tamanho da Fonte dos Labels (2025-11-25)

### Mudança

Como os nomes completos das ontologias são longos (ex: "hair follicle dermal papilla cell"), o tamanho da fonte dos labels do eixo Y foi reduzido pela metade.

### Antes e Depois

| Função | Fontsize Anterior | Fontsize Atual |
|--------|-------------------|----------------|
| `plot_comparison()` | 9 | 5 |
| `plot_individual_comparison()` | 10 | 5 |
| `plot_raw_data()` | 9 | 5 |

### Razão

Os nomes completos dos biosamples são mais descritivos mas também mais longos:
- "melanocyte of skin" 
- "hair follicle dermal papilla cell"
- "hair follicular keratinocyte"

Com fontsize menor (5 em vez de 9-10), os labels ficam mais legíveis e não se sobrepõem aos gráficos.

### Localização das Mudanças

```python
# Em todas as 3 funções de plot:
ax.set_ylabel(ylabel, fontsize=5, fontweight='bold')  # Antes: fontsize=9 ou 10
```

