# Ontology Terms Dinâmicos dos Metadados do Dataset

## Mudança Implementada

Removido o campo `ontology_terms` dos arquivos YAML de configuração. Agora as ontologias são lidas dinamicamente dos metadados do dataset.

## Problema Anterior

**Antes**: Ontologias estavam hardcoded nos arquivos YAML:

```yaml
alphagenome_api:
  enabled: true
  api_key: null
  rate_limit_delay: 0.5
  ontology_terms: ["CL:1000458", "CL:0000346", "CL:2000092"]  # ❌ Hardcoded
```

**Limitações**:
- Ontologias precisavam ser especificadas manualmente em cada config
- Risco de inconsistência entre config e dataset real
- Duplicação de informação (já está nos metadados do dataset)

## Solução Implementada

**Agora**: Ontologias lidas dinamicamente dos metadados do dataset:

```python
def get_ontology_terms_from_metadata(dataset_metadata: Dict) -> List[str]:
    """
    Extrai lista de ontology terms (CURIEs) dos metadados do dataset.
    
    Prioridade:
    1. ontology_details (contém informações completas)
    2. ontologies (lista simples de CURIEs)
    3. Fallback para lista padrão se nada encontrado
    """
    ontology_details = dataset_metadata.get('ontology_details', {})
    
    if ontology_details:
        return sorted(ontology_details.keys())
    
    ontologies = dataset_metadata.get('ontologies', [])
    if ontologies:
        return sorted([ont.strip() for ont in ontologies if ont.strip()])
    
    # Fallback
    console.print("[yellow]⚠ Ontologias não encontradas, usando lista padrão[/yellow]")
    return ["CL:1000458", "CL:0000346", "CL:2000092"]
```

**Nos arquivos YAML**:

```yaml
alphagenome_api:
  enabled: true
  api_key: null
  rate_limit_delay: 0.5
  # Nota: ontology_terms agora são lidos dinamicamente dos metadados do dataset
```

## Fonte de Dados

As ontologias são extraídas do campo `ontology_details` nos metadados do dataset:

```json
"ontology_details": {
    "CL:0000346": {
        "biosample_name": "hair follicle dermal papilla cell",
        ...
    },
    "CL:1000458": {
        "biosample_name": "melanocyte of skin",
        ...
    },
    "CL:2000092": {
        "biosample_name": "hair follicular keratinocyte",
        ...
    }
}
```

## Benefícios

1. **Sincronização automática**: Ontologias sempre correspondem ao dataset real
2. **Sem duplicação**: Informação mantida apenas nos metadados do dataset
3. **Manutenibilidade**: Não precisa editar configs ao mudar ontologias
4. **Consistência garantida**: Impossível ter mismatch entre config e dataset
5. **Ordem alfabética**: Ontologias sempre ordenadas de forma consistente

## Arquivos Modificados

### Código Python
- `neural_ancestry_predictor/verify_processed_dataset.py`
  - Adicionada função `get_ontology_terms_from_metadata()`
  - Modificadas funções: `_load_from_alphagenome_api()`, `load_raw_alphagenome_data()`
  - Todas as referências a `config['alphagenome_api']['ontology_terms']` substituídas

### Arquivos YAML (campo removido)
- `configs/verify_tyr_only.yaml`
- `configs/verify_raw_test.yaml`
- `configs/verify_processed_dataset.yaml`
- `configs/verify_alphagenome_comparison.yaml`
- `configs/verify_api_test.yaml`

## Compatibilidade

✅ **Retrocompatível**: Fallback para lista padrão se metadados não tiverem ontologias
✅ **Sem quebras**: Funciona com datasets antigos e novos
✅ **Validação**: Exibe warning se usar fallback

## Teste de Verificação

```bash
✅ Ontology terms extraídos de ontology_details:
   ['CL:0000346', 'CL:1000458', 'CL:2000092']

Hardcoded antigo (do YAML): ['CL:1000458', 'CL:0000346', 'CL:2000092']
Novo (dos metadados):        ['CL:0000346', 'CL:1000458', 'CL:2000092']

✅ Mesmas ontologias! Substituição funcionou corretamente.
```

**Nota**: A ordem mudou para alfabética (ordenada pelas chaves do dicionário).

## Integração Completa

Pipeline totalmente integrado:

1. `build_non_longevous_dataset/dataset_builder.py` → **Gera** `ontology_details` ✅
2. `neural_ancestry_predictor/verify_processed_dataset.py` → **Lê** `ontology_details` ✅
3. Configs YAML → Não precisam mais especificar ontologias ✅

## Como Funciona

1. Ao executar `verify_processed_dataset.py`, o programa carrega `dataset_metadata.json`
2. A função `get_ontology_terms_from_metadata()` extrai as ontologias dos metadados
3. As ontologias são passadas para as funções do AlphaGenome API
4. Nenhuma configuração manual necessária nos arquivos YAML

## Requer Metadados Atualizados

Para aproveitar essa funcionalidade, o dataset deve ter sido gerado com a versão atualizada do `dataset_builder.py` que inclui o campo `ontology_details`.

Para datasets antigos, a função usa fallback com lista padrão e exibe um warning.
