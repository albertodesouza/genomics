# Constantes Globais para Tamanhos de Janela

## Mudança Implementada

Substituídas todas as ocorrências hardcoded de tamanhos de janela (como `524288`, `1048576`) por constantes globais bem nomeadas no arquivo `verify_processed_dataset.py`.

## Motivação

**Problema**: Valores numéricos "mágicos" como `524288` e `1048576` apareciam repetidos em múltiplos lugares do código, tornando difícil:
- Entender o significado dos números
- Manter consistência entre diferentes partes do código
- Atualizar valores quando necessário

**Solução**: Criar constantes globais com nomes descritivos que são usadas em todo o código.

## Constantes Definidas

```python
# Tamanhos de janela suportados pelo AlphaGenome (em pares de base)
WINDOW_SIZE_32KB = 32768
WINDOW_SIZE_64KB = 65536
WINDOW_SIZE_128KB = 131072
WINDOW_SIZE_256KB = 262144
WINDOW_SIZE_512KB = 524288   # AlphaGenome: "500KB" = 524288 bp (512 KB real)
WINDOW_SIZE_1MB = 1048576

# Mapeamento de nomes de sequência para tamanhos (usado para validação)
SEQUENCE_LENGTH_MAP = {
    "SEQUENCE_LENGTH_32KB": WINDOW_SIZE_32KB,
    "SEQUENCE_LENGTH_64KB": WINDOW_SIZE_64KB,
    "SEQUENCE_LENGTH_128KB": WINDOW_SIZE_128KB,
    "SEQUENCE_LENGTH_256KB": WINDOW_SIZE_256KB,
    "SEQUENCE_LENGTH_500KB": WINDOW_SIZE_512KB,
    "SEQUENCE_LENGTH_512KB": WINDOW_SIZE_512KB,  # Alias
    "SEQUENCE_LENGTH_1MB": WINDOW_SIZE_1MB,
}

# Mapeamento reverso: tamanho -> nome de sequência
SIZE_TO_SEQUENCE_NAME = {
    WINDOW_SIZE_32KB: "SEQUENCE_LENGTH_32KB",
    WINDOW_SIZE_64KB: "SEQUENCE_LENGTH_64KB",
    WINDOW_SIZE_128KB: "SEQUENCE_LENGTH_128KB",
    WINDOW_SIZE_256KB: "SEQUENCE_LENGTH_256KB",
    WINDOW_SIZE_512KB: "SEQUENCE_LENGTH_500KB",
    WINDOW_SIZE_1MB: "SEQUENCE_LENGTH_1MB",
}
```

## Locais Modificados

### Antes (Hardcoded)
```python
window_size_map = {
    "SEQUENCE_LENGTH_500KB": 524288,
    "SEQUENCE_LENGTH_1MB": 1048576,
}
```

### Depois (Usando Constantes)
```python
window_size_map = {
    "SEQUENCE_LENGTH_500KB": WINDOW_SIZE_512KB,
    "SEQUENCE_LENGTH_1MB": WINDOW_SIZE_1MB,
}
```

## Benefícios

1. **Legibilidade**: `WINDOW_SIZE_512KB` é mais claro que `524288`
2. **Manutenibilidade**: Alterar um tamanho requer mudança em apenas um lugar
3. **Consistência**: Garante que o mesmo tamanho seja usado de forma idêntica em todo o código
4. **Documentação**: Os nomes das constantes servem como documentação inline
5. **Prevenção de erros**: Menos chances de typos ou valores incorretos

## Uso dos Mapeamentos

### SEQUENCE_LENGTH_MAP
Usado quando você tem uma string como `"SEQUENCE_LENGTH_500KB"` e precisa do tamanho em bp:

```python
target_length = SEQUENCE_LENGTH_MAP.get(window_size_key)
# Se window_size_key = "SEQUENCE_LENGTH_500KB"
# Então target_length = 524288
```

### SIZE_TO_SEQUENCE_NAME
Usado quando você tem um tamanho em bp e precisa do nome da sequência:

```python
window_size_key = SIZE_TO_SEQUENCE_NAME.get(center_bp)
# Se center_bp = 524288
# Então window_size_key = "SEQUENCE_LENGTH_500KB"
```

## Locais Substituídos

Total de substituições realizadas:
- **8 ocorrências** de `524288` substituídas por `WINDOW_SIZE_512KB`
- **11 ocorrências** de `1048576` substituídas por `WINDOW_SIZE_1MB`
- **Múltiplos dicionários** simplificados usando `SEQUENCE_LENGTH_MAP`

## Teste de Verificação

```bash
✅ Constantes definidas:
  WINDOW_SIZE_512KB = 524,288
  WINDOW_SIZE_1MB = 1,048,576

✅ SEQUENCE_LENGTH_MAP:
  SEQUENCE_LENGTH_500KB: 524,288
  SEQUENCE_LENGTH_1MB: 1,048,576

✅ SIZE_TO_SEQUENCE_NAME:
  524,288: SEQUENCE_LENGTH_500KB
  1,048,576: SEQUENCE_LENGTH_1MB

✅ Todas as constantes funcionam corretamente!
```

## Arquivo Modificado

- `neural_ancestry_predictor/verify_processed_dataset.py`

## Compatibilidade

✅ **Totalmente compatível**: Valores numéricos permanecem idênticos
✅ **Sem quebras**: Comportamento do código é 100% preservado
✅ **Melhor qualidade**: Código mais limpo e profissional

## Nota sobre AlphaGenome

AlphaGenome usa nomenclatura não convencional:
- "500KB" = 524,288 bp (que é 512 KB real em potências de 2)
- "1MB" = 1,048,576 bp (que é 1 MB exato em potências de 2)

As constantes refletem os tamanhos REAIS em bytes, não os rótulos do AlphaGenome.
