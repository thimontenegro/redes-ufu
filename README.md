# Pipeline de Rede de Proteinas (6B1T)

Implementacao em Python para:
- leitura estrutural (`.pdb`/`.cif`);
- construcao de grafo de residuos;
- analise topologica e deteccao de comunidades;
- validacao biologica (`hexon/penton`) e validacao geometrica de assembly.

## Estrutura principal

- `redes/protein_network/io.py`: parsing PDB/mmCIF e tabela de residuos
- `redes/protein_network/graph_builder.py`: construcao do grafo
- `redes/protein_network/metrics.py`: distribuicao de graus e centralidades
- `redes/protein_network/communities.py`: Louvain e Kernighan-Lin
- `redes/protein_network/validation.py`: composicao hexon/penton por comunidade
- `redes/run_pipeline.py`: pipeline principal (gera `*_summary.json` e figuras base)
- `redes/visualize_6b1t.py`: graficos/tabelas de centralidade e comunidades
- `redes/penton_assembly_validation.py`: validacao geometrica no assembly completo
- `redes/run_full_pipeline.py`: orquestra todo o fluxo

## Como executar o codigo

### 1) Execucao rapida (recomendada)

Na raiz do projeto (`Projetos`):

```bash
python redes/run_full_pipeline.py
```

O script pergunta parametros no terminal. Se apenas pressionar Enter, usa os valores padrao.

### 2) Execucao sem interacao (defaults)

```bash
python redes/run_full_pipeline.py < /dev/null
```

### 3) Execucao por etapas (opcional)

```bash
python redes/smoke_test.py
python redes/run_pipeline.py --structure redes/data/6B1T.cif --output redes/output/6b1t --cutoff 8.0 --chain-mode any --weight-mode hybrid_affinity --weighted-metrics
PYTHONPATH=redes python redes/visualize_6b1t.py --summary redes/output/6b1t/6B1T_summary.json --structure redes/data/6B1T.cif --output redes/output/6b1t/visuals
python redes/penton_assembly_validation.py --cif redes/data/6B1T.cif --output redes/output/6b1t/assembly_validation
```

## Atualizar bibliotecas/dependencias

Este projeto usa `pyproject.toml` (e `uv.lock`).

### Opcao A: usando `uv` (recomendado)

Sincronizar ambiente exatamente com o lock:

```bash
uv sync
```

Atualizar lock/dependencias para versoes mais novas:

```bash
uv lock --upgrade
uv sync
```

Adicionar nova dependencia:

```bash
uv add nome-da-biblioteca
```

### Opcao B: usando `pip`

Se estiver usando o `.venv` local:

```bash
python -m pip install --upgrade pip
python -m pip install -e .
```

Atualizar pacote especifico:

```bash
python -m pip install --upgrade pandas
```

> Observacao: se usar `pip` diretamente, o `uv.lock` pode ficar defasado do ambiente real.

## Parametros importantes

- `chain_mode`: define quais cadeias podem formar aresta.
  - `any`: permite arestas intra e intercadeia (recomendado para capsideo completo).
  - `same`: apenas mesma cadeia (analise intrachain).
  - `different`: apenas cadeias diferentes (foco em interfaces).

- `weight_mode`: define como o peso da aresta e calculado.
  - `unit`: peso fixo `1` (baseline topologico).
  - `inverse_distance`: peso `1/d` (residuos mais proximos pesam mais).
  - `hybrid_affinity`: peso `(1 + afinidade)/d` (distancia + heuristica quimica).

- `min_affinity`: afinidade minima para aceitar aresta (`0` = sem filtro quimico forte).
  - valores maiores tornam o grafo mais seletivo quimicamente.

- `weighted_metrics`: se ativado, centralidades suportadas usam o atributo `weight`.
  - em geral, melhora a leitura fisica dos caminhos quando `weight_mode != unit`.

- `connect_backbone`: conecta residuos sequenciais da mesma cadeia.
  - `True`: ajuda a evitar fragmentacao artificial da rede.
  - `False` (`--no-backbone`): usa apenas criterio de contato modelado.

### Configuracao sugerida para 6B1T

- `chain_mode=any`
- `weight_mode=hybrid_affinity`
- `min_affinity=0`
- `weighted_metrics=True`
- `connect_backbone=True`

## Saidas esperadas

Em `redes/output/6b1t`:

- `6B1T_summary.json`: resultado principal (rede, comunidades, validacao)
- `6B1T_degree_distribution.png`
- `6B1T_louvain_xy.png`
- `visuals/*`: graficos de centralidade e composicao hexon/penton por comunidade
- `assembly_validation/*`: validacao 3D de pentons/hexons no assembly completo
