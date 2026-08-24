# 02 — Briefing das Fases

> **Uma seção por fase.** O chat-ponte entrega ao chat-filho apenas o bloco da
> fase correspondente, junto com `01_Briefing_Mestre.md` e a seção "Estado atual"
> de `03_Estado_Projeto.md`.

---

## Grafo de dependências

O encadeamento entre módulos define a ordem. Nem tudo é sequencial:

```
F0 (bootstrap)
 ├─→ F1 (geodesy) ──┐
 └─→ F2 (index) ────┼─→ F3 (resolutions + cells) ─→ F4 (boundary) ─┐
      └─→ F5 (hierarchy) ─────────────────────────────────────────┼─→ F6 (topology)
                                                                   ├─→ F7 (geometry)
      └─────────────────────────────────────────────────────────→ F8 (serialization)
                                                                          │
                                                   F9 (interop + engine) ←┘
                                                             │
                                                   F10 (release v1.0.0)
```

**Paralelizáveis:** F1 e F2 após F0. F5 pode correr junto de F3 ou F4. F8 só
precisa de F2 e F5, então pode adiantar enquanto F6 e F7 correm.

**Caminho crítico:** F0 → F2 → F3 → F4 → F6 → F9 → F10.

## Quadro-resumo

| Fase | Escopo | Depende de | Risco | Status |
|---|---|---|---|---|
| F0 | Bootstrap, CI verde, `constants` + `exceptions` | — | Baixo | ⬜ |
| F1 | `geodesy` — projeção e Vincenty | F0 | Médio | ⬜ |
| F2 | `index` — parser, compose, normalize | F0 | Médio | ⬜ |
| F3 | `resolutions` + `cells` — quantização | F1, F2 | Alto | ⬜ |
| F4 | `boundary` — fronteiras e extensões | F3 | **Muito alto** | ⬜ |
| F5 | `hierarchy` — navegação hierárquica | F2 | Baixo | ⬜ |
| F6 | `topology` — vizinhança e deflexão | F4, F5 | **Muito alto** | ⬜ |
| F7 | `geometry` — polyfill, densificação | F3, F4 | Médio | ⬜ |
| F8 | `serialization` — TreeBlob + GeometryBlob | F2, F5 | Médio | ⬜ |
| F9 | `interop` + `engine` + conformidade | F6, F7, F8 | Médio | ⬜ |
| F10 | Release v1.0.0 | F9 | Baixo | ⬜ |

Legenda: ⬜ não iniciada · 🔄 em execução · ✅ concluída · ⚠️ bloqueada

---

## F0 — Bootstrap

**Objetivo.** Repositório funcional com CI verde e as constantes da especificação
transcritas.

**Escopo.**

- Merge do esqueleto no branch principal do `ICartCWB/itacart`
- `constants.py` — WGS84, tabela de resoluções, alfabetos de refinamento, zonas
  de extensão, sintaxe do índice
- `exceptions.py` — hierarquia completa de erros
- Verificar os três workflows do GitHub Actions rodando
- Configurar GitHub Pages e trusted publishing no PyPI

**Entregáveis.**

- `src/itacart/constants.py` e `exceptions.py` implementados
- CI verde nos três workflows
- `docs/project/` versionado

**Critérios de aceite.**

1. `pip install -e ".[dev]"` funciona em Linux, Windows e macOS
2. `black --check`, `isort --check-only`, `flake8`, `mypy src` limpos
3. `pytest` passa (stubs em `xfail`)
4. Toda constante da Tabela 1 do artigo confere com o PDF, valor a valor
5. `QUINARY_CODES` gera exatamente `A1`..`E5`, 25 códigos, ordem row-major

**Atenção.** A Tabela 1 tem 14 linhas; conferir uma a uma. `B-0.1` (inversão
par/ímpar) nasce de leitura apressada dessa tabela.

**Referência.** Artigo seções 3.1, 3.2; Tabela 1; Quadro 2.

---

## F1 — Geodésia

**Objetivo.** Projeção senoidal elipsoidal e geodésicas de Vincenty, sem `pyproj`
em runtime.

**Escopo.**

- `geodetic_to_sinusoidal` / `sinusoidal_to_geodetic` — Eq. (1) e (2)
- `meridian_arc` — integral elíptica por quadratura de Simpson composta,
  tolerância 1e-10 m
- `inverse_meridian_arc` — Newton sobre `meridian_arc`, derivada é o raio de
  curvatura meridiano em forma fechada
- `inverse_geodesic` / `direct_geodesic` — Vincenty 1975

**Origem.** `itacart_core/geodesy.py` (F1 + F1+ do `itacart-app`), já validado.
Portar, não reescrever.

**Entregáveis.**

- `src/itacart/geodesy.py` implementado
- `tests/unit/test_geodesy.py` com cobertura ≥ 90%
- Testes marcados `@pytest.mark.crosscheck` comparando contra `pyproj`

**Critérios de aceite.**

1. Round-trip `geodetic → sinusoidal → geodetic` fecha em < 1e-9 grau
2. `meridian_arc(90°)` confere com o quadrante meridiano do WGS84
   (10 001 965,729 m) em < 1e-4 m
3. `inverse_meridian_arc(meridian_arc(φ)) == φ` em < 1e-10 grau, para φ
   varrendo −90 a 90 em passo de 1°
4. Vincenty concorda com `pyproj.Geod` em < 1 mm para pares até 5 500 km
5. Convergência de Vincenty tratada: pontos quase antipodais levantam
   `ConvergenceError` em vez de retornar lixo

**Atenção.** A integral da Eq. (2) não tem primitiva elementar. Simpson composto
com refinamento adaptativo é o que o `itacart_core` usa e já está validado.

**Referência.** Artigo seção 2.3, Eq. (1) e (2); Zhou et al. (2007); Snyder (1987).

---

## F2 — Índice composicional

**Objetivo.** Parser, decomposição, composição e forma canônica do índice.

**Escopo.**

- `parse` — string → árvore estrutural
- `is_valid_index`, `is_atomic`, `count_cells`, `iter_cells`
- `decompose` / `compose` — a ponte entre lista plana e índice composicional
- `normalize` — forma canônica (colapso de completude + ordenação de irmãos)
- `split_components`, `quadrant_of`, `base_cell_of`

**Origem.** `itacart_core/compositional_index.py` + parser do notebook
(`DGGSTree._create_nodes`).

**Entregáveis.**

- `src/itacart/index.py` implementado
- `tests/unit/test_index.py` com cobertura ≥ 90%
- Teste de regressão para `B-0.1`

**Critérios de aceite.**

1. O índice da Figura 7 do artigo (fixture `central_park_index`) faz round-trip
   `parse → decompose → compose` sem perda
2. `normalize("...4(1,2,3,4)") == normalize("...4")` — colapso de completude
3. `normalize` é idempotente: `normalize(normalize(x)) == normalize(x)`
4. Duas escritas denotam a mesma região **se e somente se** suas formas
   canônicas são iguais — a propriedade que sustenta o Req 13 do OGC
5. `B-0.1` corrigido: resolução par aceita `1`–`4`, ímpar aceita `A1`–`E5`.
   Código fora do alfabeto do nível levanta `InvalidRefinementCodeError`
6. `decompose` tem ordem determinística (profundidade primeiro, esquerda para
   direita) — toda a semântica vetorizada depende disso

**Atenção.** `normalize` é o coração do Req 13. Sem ele o pacote não pode
reivindicar endereço único. Não confundir com `compact_cells` (F5): normalização
corrige grafia, compactação muda quais resoluções aparecem.

**Referência.** Artigo seção 3.1; exemplo `SE(1400/0374(3(C2(3))))`.

---

## F3 — Resoluções e células

**Objetivo.** Quantização (posição → célula) e geometria inversa (célula →
posição). O núcleo do DGGS.

**Escopo.**

- `resolutions.py` completo — tabela, escalas cartográficas, `refinement_ratio`
- `cells.py` — `geo_to_cell`, `sinusoidal_to_cell`, `cell_to_anchor`,
  `cell_to_centroid`, `cell_to_boundary`, `cell_to_polygon`, `cell_to_sinusoidal`
- `effective_cell_area` pode ficar como stub até F4 (depende de `cell_shape`)

**Origem.** `itacart_core/cells.py` + `sinusoidal_coordinates_to_dggs` do
notebook.

**Entregáveis.**

- `src/itacart/resolutions.py` e `cells.py` implementados
- Testes de regressão para `B-0.2`, `B-0.3`, `B-0.4`

**Critérios de aceite.**

1. Round-trip `geo_to_cell → cell_to_anchor` cai dentro da célula de origem, em
   todas as 14 resoluções
2. `cell_to_boundary` devolve 4 vértices em ordem anti-horária para célula
   interior
3. **`B-0.3` corrigido:** espelhamento correto nos quatro quadrantes. Testar
   posições simétricas — se `(10°N, 20°E)` cai em `NE(x/y)`, então `(10°S, 20°E)`
   cai em `SE(x/y)` com o mesmo par
4. **`B-0.2` corrigido:** descida em travessia única, sem recursão O(n²).
   Benchmark: resolução 13 em < 1 ms por ponto
5. **`B-0.4` resolvido:** `cell_to_centroid` existe e devolve geodésicas
6. `nominal_cell_area(r)` confere com a Tabela 1 nas 13 resoluções
7. `scale_for_resolution` confere com ambas as colunas de escala da Tabela 1
8. Área medida de célula interior (polígono do `cell_to_boundary`, reprojetado)
   bate com a nominal em < 0,01% — a propriedade de área igual

**Atenção.** A transformação paralelogramo→quadrado via residual é o truque
central; o notebook acerta isso. O espelhamento de quadrante é onde ele erra.

**Referência.** Artigo seção 3, Figura 1; Tabela 1.

---

## F4 — Comportamento de fronteira

**Objetivo.** Tornar o domínio global, completo e único. Requisitos 8–10 do OGC
Core dependem inteiramente desta fase.

**Escopo.**

- `cell_shape` — paralelogramo, triângulo, trapézio
- `is_valid_cell` — inclui a regra `X = 0` inexistente nos quadrantes ocidentais
- `is_boundary_cell`, `is_triangular_cell`, `is_trapezoidal_cell`,
  `is_equal_area_cell`
- `is_extension_cell`, `extension_zone`, `extension_zone_for_point`,
  `extension_bounds`
- `crosses_antemeridian`
- `effective_cell_area` (completar o stub de F3)
- Ajustar `geo_to_cell` e `cell_to_boundary` para os casos de fronteira

**Origem.** **Código novo.** Sem prior art no notebook nem no `itacart_core`.

**Entregáveis.**

- `src/itacart/boundary.py` implementado
- `tests/unit/test_boundary.py` com cobertura ≥ 90%
- Fixtures geográficos: Fiji, Wrangel, um ponto no meridiano de Greenwich

**Critérios de aceite.**

1. Célula que intersecta o meridiano principal em resolução 1 é triângulo
   isósceles com base = 2 × altura, espelhado em relação ao meridiano
2. Resoluções 2 a 13 **não criam** células ocidentais separadas no meridiano
   principal — são subdivisão hierárquica das células orientais adjacentes
3. `is_valid_cell` retorna `False` para resolução 1 com `X = 0` nos quadrantes
   ocidentais
4. `extension_zone_for_point` classifica corretamente: Suva (Fiji) → `"FIJI"`,
   Ilha Wrangel → `"CHUKOTKA"`, um ponto no Pacífico central → `None`
5. Limites conferem com a Figura 5: Fiji até 178°W entre 15,5°S e 21,5°S;
   Chukotka até 169,5°W entre 64°N e 72°N
6. Célula trapezoidal: `effective_cell_area` < `nominal_cell_area` e
   `is_equal_area_cell` retorna `False`
7. Toda posição em terra emergida resolve para exatamente uma célula —
   varredura sobre um conjunto de pontos de teste global

**Atenção.** Esta é a fase de maior risco do projeto. Três sub-problemas
independentes, todos código novo:

- O triângulo do meridiano principal muda a **contagem de vértices** e afeta
  `cell_to_boundary`
- As extensões mudam a **validade do domínio** e afetam `geo_to_cell`
- Os trapézios quebram a **garantia de área** e afetam `effective_cell_area`

Se a fase estourar contexto, dividir em F4a (meridiano principal) e F4b
(antemeridiano + trapézios). São separáveis.

**Referência.** Artigo seção 3.2, Figuras 4, 5 e 6.

---

## F5 — Hierarquia

**Objetivo.** Navegação hierárquica. Metade do requisito 17 do OGC Core.

**Escopo.**

- `get_parent`, `get_ancestors`, `common_ancestor`
- `get_children`, `get_descendants`, `child_position`
- `is_ancestor`, `contains`
- `compact_cells`, `uncompact_cells`

**Origem.** Análogos ASCII de `binary_index.py`; a lógica está lá em binário.

**Entregáveis.**

- `src/itacart/hierarchy.py` implementado
- `tests/unit/test_hierarchy.py` com cobertura ≥ 90%

**Critérios de aceite.**

1. `get_parent` é truncamento léxico puro — sem aritmética de ponto flutuante,
   sem índice negativo
2. `get_children` produz 4 filhos descendo para resolução par, 25 para ímpar
3. `get_descendants` é gerador, nunca lista — de resolução 1 a 13 uma célula base
   tem 10¹² descendentes
4. `compact_cells` roda até ponto fixo: célula base totalmente coberta colapsa
   até resolução 1
5. `uncompact_cells(compact_cells(x), r) == x` para `x` uniforme em resolução `r`
6. Semântica vetorizada: entrada com N células devolve N resultados alinhados
   com `decompose()`. `get_children` com `flatten=False` devolve N listas
7. `contains(region, cell)` cobre o caso em que `cell` é descendente de uma
   célula terminal da região, não só igualdade

**Atenção.** Preservar duplicatas em `get_parent` quando duas células irmãs têm o
mesmo pai — o alinhamento posicional depende disso. Quem quer deduplicar compõe
o resultado e normaliza.

**Referência.** Artigo seção 3.1.

---

## F6 — Topologia

**Objetivo.** Vizinhança e adjacência, resolvidas apenas dos índices. Completa o
requisito 17.

**Escopo.**

- `grid_disk`, `grid_ring`, `grid_distance`
- `get_neighbor`, `are_neighbor_cells`
- `deflect` — a regra de deflexão nas fronteiras
- `cells_to_directed_edge`, `directed_edge_to_cells`, `cell_to_edges`

**Origem.** **Código novo.** `B-0.5` — nada implementado em nenhuma fonte.

**Entregáveis.**

- `src/itacart/topology.py` implementado
- `tests/unit/test_topology.py` com cobertura ≥ 90%

**Critérios de aceite.**

1. Resoluções 0 e 1: aritmética inteira sobre o par `XXXX/YYYY`. Célula acima
   decrementa Y, abaixo incrementa, esquerda decrementa X, direita incrementa
2. Resoluções pares: dentro do mesmo pai, adjacência vertical soma ou subtrai 2,
   horizontal soma ou subtrai 1, sobre a numeração 1–4
3. Resoluções ímpares: dentro do mesmo pai (grade 5×5), vizinho horizontal
   incrementa o componente numérico com wrap-around entre colunas 5 e 1; vertical
   incrementa o alfabético com wrap-around entre linhas E e A
4. Vizinho em pai diferente: subir ao pai, achar o vizinho do pai, descer ao
   filho correspondente
5. `deflect` trata os três casos de fronteira: cruzar o meridiano principal cai
   em célula triangular; cruzar eixo de quadrante espelha a coordenada; cruzar o
   antemeridiano só resolve dentro de zona de extensão, senão devolve `None`
6. `are_neighbor_cells` é adjacência **por aresta** — células que se tocam apenas
   em um vértice não são vizinhas
7. `grid_disk(cell, 1, metric="chebyshev")` devolve 9 células no interior
   (origem + 8); `manhattan` devolve 5 (origem + 4)
8. Simetria: se `b ∈ grid_disk(a, 1)` então `a ∈ grid_disk(b, 1)`, inclusive
   através de fronteiras
9. Nenhuma operação chama função de coordenada — a topologia sai do índice

**Atenção.** Segunda fase de maior risco. `deflect` é a parte sutil e está
isolada de propósito, para ser testável sozinha. Documentar com clareza que
distância de grade **não é** distância geodésica — quem vem do H3 assume o
contrário.

**Referência.** Artigo seção 3.1, itens a, b, c.

---

## F7 — Geometria vetorial

**Objetivo.** Rasterização de feições, densificação e canonização.

**Escopo.**

- `polyfill` com três modos de contenção
- `count_internal_cells` — caminho rápido, acumula durante a descida
- `vertex_to_cell`, `cells_to_geometry`
- `densify_orthodromic`, `densify_segment`
- `canonicalize_rings`

**Origem.** `itacart_core/cell_filling.py`, `densification.py`,
`geometry_blob.py` (canonização) e `cadastral_processor/vertex_extractor.py`
(dedup consecutivo).

**Entregáveis.**

- `src/itacart/geometry.py` implementado
- `tests/unit/test_geometry.py` com cobertura ≥ 90%

**Critérios de aceite.**

1. `polyfill` reproduz a Figura 7 do artigo — o índice do fixture
   `central_park_index` em resoluções 6 e 7
2. Três modos de contenção se comportam como esperado:
   `contains ⊆ center ⊆ intersects`
3. `count_internal_cells` × `nominal_cell_area` bate com a área geodésica do
   polígono em < 0,1%, para polígono sem célula de fronteira
4. `count_internal_cells` não materializa índices — memória constante em
   resolução 13
5. `densify_orthodromic` produz segmentos ≤ `max_segment_m`, com pontos
   equiespaçados em distância geodésica
6. `vertex_to_cell` preserva ordem topológica; anéis mantêm o sentido de giro
7. Dedup consecutivo colapsa vértices na mesma célula; repetição **não
   consecutiva** é preservada (patologia genuína, o chamador precisa ver)
8. `canonicalize_rings` é idempotente e leva o mesmo anel iniciado em vértices
   diferentes à mesma sequência. LINESTRING é direcional e **não** rotaciona
9. Geometria cruzando o antemeridiano fora de zona de extensão levanta
   `AntemeridianError`

**Atenção.** A densificação existe porque reta no plano senoidal não é geodésica
no elipsoide. Sem ela, aresta longa preenche as células erradas no meio do
caminho.

**Referência.** Artigo seção 4, Figura 7; Tong et al. (2013).

---

## F8 — Serialização

**Objetivo.** Portar os dois formatos binários.

**Escopo.**

- `tree_blob.py` — TreeBlob completo, forma canônica prefix-index
- `geometry_blob.py` — GeometryBlob v2, sequência ordenada com header de 9 bytes

**Origem.** `itacart_core/binary_index.py` (1 217 linhas) e `geometry_blob.py`
(1 675 linhas). **Portar, não reescrever** — são 180 testes já validados.

**Entregáveis.**

- `src/itacart/serialization/` implementado
- `tests/unit/test_serialization.py` com cobertura ≥ 88%
- `docs/binary_encoding_spec.md` e `geometry_encoding_spec.md` migrados

**Critérios de aceite.**

1. Propriedades formais F-1 a F-8 do TreeBlob preservadas: determinismo,
   `decode∘encode = recompose`, `encode = encode∘recompose`, `recompose`
   idempotente, round-trip bit-exato, content-addressability, folha em
   resolução 13 = 10 bytes, `prefix_at_resolution` monotônica
2. Densidade do TreeBlob confere com o medido: ~5,13 b/v assintótico em 1 000
   vértices
3. GeometryBlob: header de 9 bytes conforme a spec; round-trip preserva ordem,
   topologia de anéis, modelo de aresta e densificação
4. `geometry_to_tree` é one-way — mesmo conjunto de vértices em sequências
   diferentes produz o mesmo TreeBlob. Isso é *by design*, não bug
5. Aplicar o patch `P-2.2.9` do `itacart-app` (§17, envelope de commitment) se
   ainda pendente
6. `B-F3.1` do `itacart-app` corrigido: `geometry_to_tree` concatenava ASCII com
   vírgulas incompatível com `encode_tree`

**Atenção.** O `itacart-app` tem uma pendência aberta (`P-F3.1`) exatamente em
`geometry_to_tree`. Verificar se foi resolvida lá antes de portar; se não,
resolver aqui e reportar de volta.

**Referência.** `docs/binary_encoding_spec.md` e `geometry_encoding_spec.md` v2
do `itacart-app`.

---

## F9 — Interoperabilidade, engine e conformidade

**Objetivo.** Fechar os requisitos declarativos do OGC e deixar a suite de
conformidade verde.

**Escopo.**

- `interop.py` — GeoJSON, WKT, GeoDataFrame (Req 18–19)
- `engine.py` — classe `ITACaRT`, `describe`, `crs`, `conformance` (Req 6–7)
- `tests/conformance/` — remover todos os `xfail`

**Entregáveis.**

- `src/itacart/interop.py` e `engine.py` implementados
- Suite de conformidade passando sem `xfail`
- Página `docs/source/concepts/conformance.rst` escrita

**Critérios de aceite.**

1. `describe()` cobre identidade, DOI, CRS, método de tesselação, geometria da
   célula, tabela de resoluções, razões de refinamento e tratamentos de fronteira
2. `conformance()` reproduz os Quadros 4 e 5 do artigo, requisito a requisito
3. Todo teste de conformidade passa — Core integral, EAERS com as divergências
   pinadas como intencionais
4. `test_req_22_25_not_met_by_design` verifica que **não** existe interface
   poliédrica; se um dia existir, o teste falha e a divergência é reavaliada
5. `cells_to_geojson` produz FeatureCollection válido em EPSG:4326
6. `to_geodataframe` funciona com `crs="EPSG:4326"` e com `SINUSOIDAL_PROJ`
7. Job `conformance` do CI verde

**Atenção.** O artigo marca Req 18–19 como "Met (by design)". Esta fase é o que
transforma isso em verificável.

**Referência.** Artigo seção 4, Quadros 4 e 5; Gibb (2021).

---

## F10 — Release v1.0.0

**Objetivo.** Publicar no PyPI com documentação completa.

**Escopo.**

- Páginas de conceito em `docs/source/concepts/` escritas, **com as figuras do
  artigo** (orientação do orientador: maximizar gráficos)
- README revisado com exemplos que rodam
- `CHANGELOG.md` consolidado
- Notebook de exemplo reproduzindo a Figura 7
- Bump para `1.0.0`, tag, release no GitHub

**Entregáveis.**

- Pacote publicado em `pypi.org/project/itacart/`
- Documentação no GitHub Pages
- Notebook em `examples/`

**Critérios de aceite.**

1. `pip install itacart` funciona em ambiente limpo
2. Todo exemplo do README executa sem erro
3. Sphinx compila com `-W` (warning é erro)
4. Cobertura total ≥ 88%
5. `twine check` limpo
6. DOI e citação BibTeX corretos no README e na documentação

**Atenção.** Publicação é irreversível — versão queimada no PyPI não volta.
Validar em TestPyPI antes.

---

## Fases candidatas pós-v1

Fora do escopo da v1, registradas para não se perderem:

- **F11** — Benchmarking quantitativo contra H3, S2 e rHEALPix. Insumo direto
  do Grid Lab e da seção de trabalhos futuros do artigo
- **F12** — Backend OGC API-DGGS via `pydggsapi`
- **F13** — Otimização de performance: descida vetorizada, cache de projeção
- **F14** — Codificação delta cross-cell para regime muito grande
  (`F2.3` parked no `itacart-app`, candidato a segundo artigo)
