# Carregamento Dinâmico de Genes dos Metadados (ATUALIZADO)

## Mudança Implementada

Modificado `verify_processed_dataset.py` para ler a lista de genes disponíveis dinamicamente dos metadados do dataset. **Se o campo 'genes' não for encontrado, o programa termina com erro.**

## Problema Anterior

**Antes**: Lista de genes estava hardcoded na classe `InteractiveComparisonViewer`:

```python
# Lista de genes disponíveis
self.available_genes = ["SLC24A5", "SLC45A2", "OCA2", "HERC2", "MC1R", 
                       "EDAR", "MFSD12", "DDB1", "TCHH", "TYR", "TYRP1"]
```

**Limitações**:
- Lista manual precisava ser atualizada para cada novo dataset
- Não refletia automaticamente mudanças nos genes do dataset
- Inconsistente se dataset tivesse genes diferentes

## Solução Implementada

**Agora**: Leitura dinâmica dos metadados do dataset com validação obrigatória:

```python
# Ler lista de genes disponíveis dos metadados do dataset
self.available_genes = dataset_metadata.get('genes', [])
if not self.available_genes:
    # Erro: campo 'genes' não encontrado nos metadados
    console.print("[red]❌ ERRO: Campo 'genes' não encontrado nos metadados do dataset![/red]")
    console.print("[yellow]Os metadados do dataset precisam ser regenerados com a versão atualizada do dataset_builder.py[/yellow]")
    console.print("[cyan]Execute: python3 build_non_longevous_dataset.py --config <seu_config.yaml>[/cyan]")
    console.print("[cyan]Com generate_dataset_metadata: true no config[/cyan]")
    raise ValueError(
        "Campo 'genes' não encontrado nos metadados do dataset. "
        "Regenere os metadados do dataset com a versão atualizada do dataset_builder.py"
    )
```

## Benefícios

1. **Dinâmico**: Reflete automaticamente os genes presentes no dataset
2. **Validação estrita**: Garante que metadados estejam atualizados
3. **Mensagens claras**: Indica exatamente como corrigir o problema
4. **Manutenção**: Não precisa editar código ao adicionar novos genes
5. **Consistência**: Lista sempre sincronizada com o dataset real
6. **Ordem alfabética**: Genes aparecem ordenados alfabeticamente

## Comportamento

### ✅ Sucesso (com campo 'genes')
```
✅ Genes encontrados: 11
   Lista: ['DDB1', 'EDAR', 'HERC2', 'MC1R', 'MFSD12', 'OCA2', 
           'SLC24A5', 'SLC45A2', 'TCHH', 'TYR', 'TYRP1']
```

### ❌ Erro (sem campo 'genes')
```
❌ ERRO: Campo 'genes' não encontrado nos metadados do dataset!
Os metadados do dataset precisam ser regenerados com a versão atualizada do dataset_builder.py
Execute: python3 build_non_longevous_dataset.py --config <seu_config.yaml>
Com generate_dataset_metadata: true no config

ValueError: Campo 'genes' não encontrado nos metadados do dataset. 
Regenere os metadados do dataset com a versão atualizada do dataset_builder.py
```

## Como Corrigir Datasets Antigos

Se você encontrar o erro acima, regenere os metadados do dataset:

### Opção 1: Regenerar apenas metadados (rápido)

```bash
cd /home/lume2/genomics/build_non_longevous_dataset
python3 -c "
from pathlib import Path
from dataset_builder import DatasetMetadataBuilder

dataset_dir = Path('/dados/GENOMICS_DATA/top3/non_longevous_results_genes')
window_size = 524288  # Ajuste conforme seu config

builder = DatasetMetadataBuilder(
    dataset_dir=dataset_dir,
    dataset_name='non_longevous_dataset_genes',
    window_size=window_size
)

builder.scan_individuals()
builder.save_metadata()
builder.print_summary()
"
```

### Opção 2: Re-executar pipeline completo

```bash
cd /home/lume2/genomics/build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/small_gene.yaml
```

Com no config:
```yaml
pipeline:
  steps:
    generate_dataset_metadata: true
```

## Dependências

**REQUER** que o dataset tenha sido gerado com a versão atualizada do `dataset_builder.py` que inclui o campo `genes` nos metadados.

**Não há mais fallback** - o campo é obrigatório.

## Arquivo Modificado

- `neural_ancestry_predictor/verify_processed_dataset.py` (linhas 957-971)

## Compatibilidade

⚠️ **NÃO retrocompatível por design**: Força atualização dos metadados
✅ **Validação estrita**: Garante consistência do sistema
✅ **Mensagens de erro claras**: Indica exatamente como resolver

## Vantagens da Validação Estrita

1. **Prevenção de inconsistências**: Não permite usar datasets com metadados desatualizados
2. **Força boas práticas**: Incentiva manter metadados atualizados
3. **Debugging mais fácil**: Erro claro em vez de comportamento silencioso incorreto
4. **Documentação através do erro**: Mensagem explica como resolver
