# Genomics Workbench: Planejamento de Ferramentas Graficas

Este documento orienta a implementacao paralela de um conjunto de ferramentas graficas locais para o repositorio. O publico alvo sao pessoas com conhecimento de biologia que precisam visualizar dados, construir novos genomas, executar pipelines e comparar experimentos sem depender fortemente do terminal.

## 1. Objetivo

Criar um workbench web local, iniciado por Python, que exponha visualizacoes e acoes seguras sobre os dados genomicos ja existentes no repositorio.

O sistema deve priorizar:

- leitura sob demanda de arquivos grandes;
- baixo uso de memoria em maquinas com 16 GB RAM;
- interfaces simples com filtros, previews e botoes de execucao;
- comandos equivalentes visiveis para reprodutibilidade;
- implementacoes independentes para permitir trabalho paralelo de IAs/agentes.

## 2. Principios Tecnicos

- Nao carregar datasets completos em memoria.
- Usar paginacao, janelas e indices auxiliares pequenos.
- Preferir servidores locais simples em Python no curto prazo.
- Evitar dependencias frontend pesadas inicialmente.
- Usar HTML/CSS/JS puro ou Plotly quando grafico interativo for necessario.
- Manter cada visualizacao em arquivo proprio ate estabilizar.
- Nao modificar pipelines centrais sem necessidade.
- Preservar compatibilidade com o ambiente `genomics` ativado por `scripts/start_genomics_universal.sh`.
- Sempre expor limites configuraveis (`--max-rows`, `--max-cells`, `--page-size`).
- Mostrar erros de forma amigavel na UI, mas manter stack trace no console quando util.

## 3. Estrutura Recomendada

Implementacao incremental em arquivos independentes:

```text
genotype_based_predictor/
  aligned_dna_viewer.py          # ja iniciado
  dataset_browser.py             # browser de dataset/metadados
  view_builder.py                # construtor de views JSON
  experiment_dashboard.py        # metricas e comparacao de experimentos
  pipeline_runner.py             # execucao de pipelines/jobs locais
  alphagenome_track_viewer.py    # visualizacao de tracks AlphaGenome
  tensor_viewer.py               # inspecao de tensores processados/cache
  model_error_viewer.py          # analise de erros do modelo
  ui_common.py                   # opcional: helpers compartilhados
  docs/
    genomics_workbench_plan.md
```

Se varios agentes trabalharem em paralelo, cada agente deve editar preferencialmente apenas seu arquivo e, quando necessario, adicionar secoes neste documento.

### 3.1 Launcher Unificado com Roteamento

O arquivo `genotype_based_predictor/genomics_workbench.py` inicia as visualizacoes implementadas como subprocessos e serve uma pagina principal com roteamento em uma unica porta. Internamente, cada viewer continua rodando como servidor isolado; externamente, o usuario acessa apenas o Workbench.

Comando recomendado:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.apps.genomics_workbench
```

Depois abra:

```text
http://127.0.0.1:8780
```

Rotas publicas:

```text
http://127.0.0.1:8780/apps/dataset/
http://127.0.0.1:8780/apps/view-builder/
http://127.0.0.1:8780/apps/experiments/
http://127.0.0.1:8780/apps/tracks/
http://127.0.0.1:8780/apps/alignment/
```

Portas internas padrao:

- Dataset Browser em `http://127.0.0.1:8770`;
- View Builder em `http://127.0.0.1:8771`;
- Experiment Dashboard em `http://127.0.0.1:8772`;
- AlphaGenome Track Viewer em `http://127.0.0.1:8774`;
- Aligned DNA Viewer em `http://127.0.0.1:8765`.

Essas portas internas sao detalhes de implementacao. A pagina principal atua como reverse proxy local para que as aplicacoes funcionem pela mesma porta (`8780`) com rotas `/apps/<app>/`. Ao navegar entre cards, a visualizacao correspondente ja esta em execucao.

Parametros uteis:

```bash
python3 -m genotype_based_predictor.apps.genomics_workbench \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --runs-root /dados/GENOMICS_DATA/v1/1kG_high_coverage_runs \
  --aligned-tsv-root genotype_based_predictor/aligned_dna_genes_1000_all \
  --host 127.0.0.1 \
  --port 8780
```

Logs dos subprocessos ficam em `/tmp/genomics_workbench_logs` por padrao, ou no diretorio informado com `--log-dir`.

## 4. Contratos Comuns

### 4.1 Execucao

Cada visualizacao deve ser executavel como modulo:

```bash
source scripts/start_genomics_universal.sh && python3 -m genotype_based_predictor.<modulo> <args>
```

Cada servidor deve aceitar:

- `--host`, default `127.0.0.1`;
- `--port`, default distinto por modulo;
- limites de seguranca quando aplicavel.

### 4.2 Respostas JSON

APIs internas devem retornar JSON com:

```json
{"error": "mensagem"}
```

em caso de erro tratavel.

### 4.3 Arquivos Grandes

Para TSVs ou arquivos linha-orientados grandes:

- criar indice `.idx.json` com offsets;
- guardar `mtime_ns` e `size` para invalidar cache;
- ler apenas as linhas ou janelas solicitadas.

Para arrays/tensores:

- carregar apenas shards ou amostras selecionadas;
- evitar empilhar batches inteiros;
- calcular estatisticas em streaming quando possivel.

### 4.4 UI

Cada UI deve ter:

- formulario de filtros;
- area de status/erro;
- tabela ou grafico principal;
- explicacao curta dos limites;
- destaque visual consistente.

Paleta sugerida:

- fundo escuro;
- vermelho para mutacao/erro;
- amarelo para gap/aviso;
- azul para links/acoes;
- verde para sucesso.

## 5. Visualizacoes Planejadas

## 5.1 Aligned DNA Viewer

Status: implementado em `aligned_dna_viewer.py` e integrado ao Workbench em `/apps/alignment/`.

Objetivo: visualizar TSVs alinhados gerados por `export_aligned_dna.py`.

Funcionalidades atuais:

- lista genes a partir de arquivos `.tsv`;
- seleciona individuos por checkboxes com filtros de superpopulacao, populacao e busca textual;
- seleciona haplotipo (`H1`, `H2`, ambos);
- escolhe janela (`start`, `length`);
- possui botoes de janela anterior/proxima, incluindo botao inferior ao fim da tabela;
- destaca mutacoes em relacao a referencia;
- destaca `X`;
- filtra apenas posicoes variantes;
- usa indice por offset.
- mostra coluna `index` com a posicao no eixo expandido;
- mostra coluna `ref genome pos` com a coordenada genomica de referencia;
- deixa `ref genome pos` vazio quando o `REF` contem `X` naquela linha;
- mostra variantes do VCF que caem na janela selecionada.

Melhorias futuras:

- exportar janela visivel como TSV/PNG;
- mini-mapa de variantes;
- links cruzados com tracks AlphaGenome.

Documentacao operacional atual:

```text
genotype_based_predictor/README.md
```

## 5.2 Dataset Browser

Arquivo sugerido: `dataset_browser.py`.

Objetivo: navegar metadados do dataset e identificar completude/problemas.

Entradas:

- `dataset_dir` contendo `dataset_metadata.json`, `individuals/`, `references/`.

Funcionalidades MVP:

- listar individuos paginados;
- filtrar por `sample_id`, `population`, `superpopulation`, `sex`;
- listar genes disponiveis;
- mostrar distribuicao por populacao/superpopulacao;
- mostrar completude individuo x gene em modo resumido;
- abrir detalhe de individuo com metadados e janelas;
- API `/api/summary`, `/api/individuals`, `/api/genes`, `/api/individual/<id>`.

Cuidados:

- nao varrer todos os arquivos profundamente por padrao;
- usar metadados globais quando disponiveis;
- completude detalhada deve ser sob demanda.

## 5.3 View Builder

Arquivo sugerido: `view_builder.py`.

Objetivo: criar arquivos `.view.json` sem editar JSON manualmente.

Entradas:

- `dataset_dir`;
- opcionalmente uma view existente para clonar.

Funcionalidades MVP:

- selecionar genes;
- selecionar individuos por filtros;
- selecionar superpopulacoes/populacoes;
- configurar `alphagenome_outputs`, `window_center_size`, `normalization_method`, `ontology_terms`;
- mostrar preview com numero de genes e amostras;
- salvar JSON em `genotype_based_predictor/views/<nome>.view.json`;
- gerar snippet YAML para `dataset_input.view_path`.

Cuidados:

- validar paths;
- nao sobrescrever view existente sem confirmacao;
- usar schema compativel com `PipelineConfig`.

## 5.4 Experiment Dashboard

Arquivo sugerido: `experiment_dashboard.py`.

Objetivo: explorar experimentos salvos em `processed_cache_dir`.

Entradas:

- diretorio base de runs, por exemplo `/dados/GENOMICS_DATA/v1/1kG_high_coverage_runs`.

Funcionalidades MVP:

- listar experimentos;
- detectar `config.yaml`, `models/*.pt`, `*_results.json`;
- extrair metricas (`accuracy`, `precision`, `recall`, `f1`);
- ordenar por metrica;
- visualizar matriz de confusao quando existir;
- comparar configs em texto/diff simples;
- API `/api/experiments`, `/api/experiment/<name>`, `/api/results/<name>`.

Cuidados:

- nao carregar checkpoints `.pt` por padrao;
- ler apenas JSON/YAML pequenos;
- tratar experimentos incompletos.

## 5.5 Pipeline Runner

Arquivo sugerido: `pipeline_runner.py`.

Objetivo: executar pipelines comuns pela UI.

Funcionalidades MVP:

- selecionar um comando predefinido;
- preencher argumentos;
- mostrar comando equivalente;
- executar subprocesso;
- streaming de logs;
- status de jobs;
- cancelar job;
- historico em JSON local.

Comandos iniciais:

- gerar TSV alinhado (`export_aligned_dna.py`);
- gerar viewer columns (`aligned_dna_columns.py`);
- iniciar treino (`train.py`);
- iniciar web viewers;
- gerar alinhamentos para todos os genes de uma view.

Observacao: a geracao de TSVs para todos os individuos deve usar `export_aligned_dna.py --batch-size` para limitar memoria. Em maquinas com 16 GB RAM, usar `--batch-size 2` ou `--batch-size 4` e evitar varios genes em paralelo.

Cuidados:

- usar allowlist de comandos;
- nao aceitar shell arbitrario no MVP;
- registrar cwd, timestamp e exit code;
- evitar multiplos jobs pesados simultaneos por padrao.

## 5.6 AlphaGenome Track Viewer

Arquivo sugerido: `alphagenome_track_viewer.py`.

Objetivo: visualizar arrays de predicao AlphaGenome por individuo/gene/haplotipo/output.

Funcionalidades MVP:

- selecionar individuo, gene, haplotipo, output;
- listar tracks/metadata quando disponivel;
- renderizar janela de sinal como linha ou heatmap;
- comparar H1 vs H2;
- mostrar estatisticas min/max/media;
- API com amostragem/downsampling para nao enviar arrays enormes.

Cuidados:

- carregar somente arquivo de predicao solicitado;
- downsample no backend para plotagens longas;
- limitar numero de tracks simultaneas.

## 5.7 Tensor Viewer

Arquivo sugerido: `tensor_viewer.py`.

Objetivo: inspecionar tensores processados/cacheados que entram no modelo.

Funcionalidades MVP:

- selecionar cache de dataset;
- selecionar split (`train`, `val`, `test`);
- selecionar indice/amostra;
- carregar tensor de um shard;
- mostrar shape, target, classe;
- heatmap com downsampling;
- estatisticas por canal;
- mostrar canais de sinal/mascaras quando aplicavel.

Cuidados:

- ler `shards_index.json`;
- carregar apenas shard necessario;
- nao carregar split inteiro.

## 5.8 Model Error Viewer

Arquivo sugerido: `model_error_viewer.py`.

Objetivo: visualizar erros de classificacao a partir dos resultados salvos.

Funcionalidades MVP:

- carregar `*_results.json`;
- mostrar matriz de confusao clicavel;
- listar classes com maior erro;
- se houver predicoes por amostra no futuro, listar amostras erradas;
- links para abrir individuo/gene nos viewers.

Cuidados:

- atualmente `evaluation.py` salva metricas agregadas, nao predicoes por amostra;
- documentar necessidade futura de salvar `predictions.csv`.

## 6. Plano de Paralelizacao Para Agentes

Cada agente deve trabalhar em um arquivo independente e retornar:

- arquivos criados/modificados;
- comando de smoke test;
- limitacoes conhecidas;
- proximos passos.

### Agente A: Dataset Browser

Escopo:

- implementar `genotype_based_predictor/dataset_browser.py`;
- UI web local;
- APIs de resumo e individuos paginados;
- leitura de `dataset_metadata.json` e metadados individuais sob demanda.

Smoke test sugerido:

```bash
python3 -m genotype_based_predictor.apps.dataset_browser /dados/GENOMICS_DATA/v1/1kG_high_coverage --port 8770
```

### Agente B: View Builder

Escopo:

- implementar `genotype_based_predictor/view_builder.py`;
- criar views JSON via UI;
- usar dataset metadata para genes e individuos;
- salvar arquivo em path informado.

Smoke test sugerido:

```bash
python3 -m genotype_based_predictor.apps.view_builder /dados/GENOMICS_DATA/v1/1kG_high_coverage --port 8771
```

### Agente C: Experiment Dashboard

Escopo:

- implementar `genotype_based_predictor/experiment_dashboard.py`;
- listar experimentos e metricas;
- visualizar matriz de confusao.

Smoke test sugerido:

```bash
python3 -m genotype_based_predictor.apps.experiment_dashboard /dados/GENOMICS_DATA/v1/1kG_high_coverage_runs --port 8772
```

### Agente D: Pipeline Runner

Escopo:

- implementar `genotype_based_predictor/pipeline_runner.py`;
- UI para jobs allowlisted;
- streaming/polling de logs;
- cancelar job.

Smoke test sugerido:

```bash
python3 -m genotype_based_predictor.pipeline_runner --port 8773
```

### Agente E: AlphaGenome Track Viewer

Escopo:

- implementar `genotype_based_predictor/alphagenome_track_viewer.py`;
- visualizar predicoes por individuo/gene/haplotipo/output;
- downsampling no backend.

Smoke test sugerido:

```bash
python3 -m genotype_based_predictor.apps.alphagenome_track_viewer /dados/GENOMICS_DATA/v1/1kG_high_coverage --port 8774
```

### Agente F: Tensor Viewer

Escopo:

- implementar `genotype_based_predictor/tensor_viewer.py`;
- abrir cache processado;
- carregar uma amostra/shard por vez;
- renderizar heatmap downsampled.

Smoke test sugerido:

```bash
python3 -m genotype_based_predictor.tensor_viewer /dados/GENOMICS_DATA/v1/1kG_high_coverage_runs/datasets --port 8775
```

## 7. Roadmap Incremental

### Fase 1: Ferramentas independentes

- Dataset Browser MVP.
- Aligned DNA Viewer MVP.
- Experiment Dashboard MVP.
- View Builder MVP.

### Fase 2: Execucao e integracao

- Pipeline Runner MVP.
- Links cruzados entre viewers.
- Historico de jobs.

### Fase 3: Dados numericos e modelo

- AlphaGenome Track Viewer.
- Tensor Viewer.
- Model Error Viewer.

### Fase 4: Workbench unificado

- Servidor unico com abas.
- Autenticacao local opcional.
- Relatorios exportaveis.
- Persistencia em SQLite.

## 8. Riscos

- Arquivos grandes podem travar o navegador se limites nao forem aplicados.
- Jobs longos podem deixar subprocessos orfaos se cancelamento nao for bem tratado.
- Paths absolutos do dataset podem variar entre maquinas.
- Alguns caches podem estar incompletos ou em formato legado.
- Resultado de modelos atualmente nao inclui predicao por amostra, limitando analise de erros.

## 9. Definicao de Pronto Para Cada Viewer

Um viewer MVP esta pronto quando:

- roda com `python3 -m genotype_based_predictor.<modulo>`;
- abre no navegador local;
- possui pelo menos uma API JSON testavel;
- nao carrega tudo em memoria por padrao;
- mostra erros amigaveis;
- possui comando de smoke test documentado;
- funciona com dados reais ou falha com mensagem clara se o dado esperado nao existe.
