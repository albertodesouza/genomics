# Melhorias no verify_processed_dataset.py (ATUALIZADO)

## Resumo das Mudanças

Modificado `verify_processed_dataset.py` para ler a lista de genes dinamicamente dos metadados do dataset. **O campo 'genes' é agora obrigatório - o programa termina com erro se não for encontrado.**

## Mudança Específica

### Localização
Classe `InteractiveComparisonViewer`, método `__init__` (linhas ~957-971)

### Antes
```python
# Lista de genes disponíveis
self.available_genes = ["SLC24A5", "SLC45A2", "OCA2", "HERC2", "MC1R", 
                       "EDAR", "MFSD12", "DDB1", "TCHH", "TYR", "TYRP1"]
```

### Depois
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

## Vantagens

✅ **Dinâmico**: Lista de genes sincronizada automaticamente com o dataset  
✅ **Manutenível**: Não requer edição de código ao adicionar genes  
✅ **Validação estrita**: Garante que metadados estejam sempre atualizados  
✅ **Mensagens claras**: Erro explica exatamente como resolver o problema  
✅ **Consistente**: Ordem alfabética dos metadados  

## Comportamento

### Com campo 'genes' presente ✅
```bash
✅ Genes encontrados: 11
   Lista: ['DDB1', 'EDAR', 'HERC2', 'MC1R', 'MFSD12', 'OCA2', 
           'SLC24A5', 'SLC45A2', 'TCHH', 'TYR', 'TYRP1']
```

### Sem campo 'genes' ❌
```bash
❌ ERRO: Campo 'genes' não encontrado nos metadados do dataset!
Os metadados do dataset precisam ser regenerados com a versão atualizada do dataset_builder.py
Execute: python3 build_non_longevous_dataset.py --config <seu_config.yaml>
Com generate_dataset_metadata: true no config

ValueError: Campo 'genes' não encontrado nos metadados do dataset.
```

## Integração com build_non_longevous_dataset

Esta mudança complementa as melhorias feitas em `build_non_longevous_dataset`:

1. **dataset_builder.py** agora gera campo `genes` nos metadados ✅
2. **verify_processed_dataset.py** agora **requer** campo `genes` nos metadados ✅

Resultado: Pipeline totalmente integrado com validação estrita.

## Compatibilidade

⚠️ **NÃO retrocompatível por design**: Força atualização dos metadados  
✅ **Validação garantida**: Previne uso de datasets com metadados desatualizados  
✅ **Erro auto-documentado**: Mensagem explica como corrigir  

## Como Corrigir Datasets Antigos

Se encontrar o erro, regenere os metadados:

```bash
cd /home/lume2/genomics/build_non_longevous_dataset
python3 -c "
from pathlib import Path
from dataset_builder import DatasetMetadataBuilder

dataset_dir = Path('/path/to/your/dataset')
window_size = 524288  # Ajuste conforme necessário

builder = DatasetMetadataBuilder(
    dataset_dir=dataset_dir,
    dataset_name='your_dataset_name',
    window_size=window_size
)

builder.scan_individuals()
builder.save_metadata()
"
```

## Documentação Completa

Veja `DYNAMIC_GENE_LOADING_UPDATED.md` para detalhes técnicos e exemplos completos.
