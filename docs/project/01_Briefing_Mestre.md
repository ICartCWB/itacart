# 01 — Briefing Mestre

> **Documento de contexto.** Todo chat-filho recebe este arquivo na abertura da
> sua fase. Descreve o que é o pacote, por que existe, e quais decisões já foram
> tomadas e não devem ser re-litigadas.

---

## 1. O que é o pacote `itacart`

Implementação de referência em Python do **ITACaRT** (ITA Cadastral Ellipsoidal
Reference Tessellation), um Discrete Global Grid System de células paralelogramas
de área igual, tesselado diretamente sobre o elipsoide WGS84, projetado para
mapeamento cadastral terrestre.

- **PyPI:** `pypi.org/project/itacart/` (placeholder em 0.1.0a3)
- **Repositório:** `github.com/ICartCWB/itacart`
- **Artigo canônico:** DOI 10.14393/rbcv77n0a-79281 — *Revista Brasileira de
  Cartografia*, v. 77, 2025
- **Autores:** Israel Nunes da Silva (ITA), Gabriel Dietzsch (IEAv),
  Elcio Hideiti Shiguemori (IEAv)

O artigo é a especificação. Quando o código e o artigo divergirem, o artigo vence.

## 2. Por que existe

DGGS consagrados foram projetados para outros fins: H3 para logística em escala
web, S2 para indexação de banco de dados, GEOSOT para alinhamento com a grade de
lat/lon. Nenhum prioriza simultaneamente **isometria de área** e **precisão
geodésica**, que são requisitos não-negociáveis do cadastro — a área da parcela é
atributo jurídico e fiscal.

O rHEALPix chega perto, mas sua resolução máxima (aresta ~2 m) não atende
levantamento cadastral.

ITACaRT preenche essa lacuna com quatro escolhas deliberadas:

1. **Tesselação direta no elipsoide**, sem poliedro intermediário — fidelidade
   geodésica acima de conveniência computacional
2. **Células paralelogramas** em vez de quadrados — absorve melhor a distorção
   angular inerente à projeção senoidal
3. **Hierarquia decimal** (14 níveis, 10 km a 1 cm) com áreas métricas inteiras —
   habilita tokenização onde 1 token = 1 unidade padrão de área
4. **Índice composicional** — representa uma feição vetorial inteira em uma
   string, em vez de uma lista desestruturada de identificadores atômicos

## 3. Proveniência do código

Três fontes alimentam este pacote. Saber de onde vem cada parte evita
reimplementar o que já foi validado.

### 3.1 `itacart_core/` do projeto `itacart-app`

Motor DGGS em Python puro, **406 testes, ~99% de cobertura**, `shapely` como
única dependência de runtime. Reproduz a Figura 7 do artigo bit-a-bit.

Módulos aproveitáveis:

| Origem | Destino no pacote | Observação |
|---|---|---|
| `geodesy.py` | `geodesy.py` | Senoidal elipsoidal + Vincenty; 12 µm vs `pyproj` |
| `resolutions.py` | `resolutions.py` | Tabela de resoluções |
| `cells.py` | `cells.py` | `point_to_cell`, `cell_to_point` |
| `compositional_index.py` | `index.py` | Parser ASCII |
| `cell_filling.py` | `geometry.py` | Descida hierárquica sobre plano senoidal |
| `densification.py` | `geometry.py` | `densify_orthodromic` |
| `binary_index.py` | `serialization/tree_blob.py` | TreeBlob, 1 217 linhas |
| `geometry_blob.py` | `serialization/geometry_blob.py` | GeometryBlob v2, 1 675 linhas |
| `engine.py` | `engine.py` | `IDGGSEngine` |

**Não migram:** `cadastral_processor/` e `blockchain_client/` são específicos da
FAB. A lógica genérica de `vertex_extractor.py` (densificação, dedup consecutivo)
migra para `geometry.py`; o resto fica fora.

### 3.2 Notebook `DGGS_Tree.ipynb`

Protótipo em Colab. Bom para o parser do índice composicional e para a
geocodificação coordenada→índice (transformação paralelogramo→quadrado via
residual). **Contém bugs conhecidos** — ver seção 5.

### 3.3 Código novo

Sem prior art em nenhuma das fontes:

- **Comportamento de fronteira** (`boundary.py`) — triângulos no meridiano
  principal, extensões do antemeridiano, trapézios
- **Vizinhança** (`topology.py`) — `grid_disk`, deflexão entre quadrantes
- **Interoperabilidade** (`interop.py`) — GeoJSON, WKT, GeoDataFrame
- **Suite de conformidade OGC** (`tests/conformance/`)

## 4. Decisões de projeto já tomadas

Não re-litigar sem motivo forte. Cada uma tem alternativa rejeitada registrada.

### D-0.1 — Escopo completo na v1

Comportamento de fronteira e vizinhança entram na v1, não em versão futura. São
parte da definição do DGGS: sem eles os requisitos 8–10 do OGC Core (domínio
global, completo e único) não são atendidos, e o pacote não pode reivindicar
conformidade.

*Rejeitado:* v1 mínima com fronteiras adiadas.

### D-0.2 — Semântica vetorizada em vez de tipos separados

Toda operação aceita índice composicional. Operações que retornam um valor por
célula retornam lista alinhada posicionalmente com `decompose()`.

*Rejeitado:* tipos `Cell` e `Region` distintos. O índice composicional já é
intrinsecamente um conjunto ordenado; bifurcar em dois tipos duplica a superfície
sem ganho.

### D-0.3 — Sem `pyproj` em runtime

A projeção senoidal elipsoidal é implementada diretamente das Eq. (1) e (2) do
artigo. `pyproj` fica como dependência de desenvolvimento, apenas para
verificação cruzada em testes.

*Rejeitado:* depender de PROJ. Dependência pesada para duas equações fechadas, e
tira o controle sobre a integral elíptica.

### D-0.4 — `cell_to_anchor` separado de `cell_to_centroid`

São funções distintas. O *anchor* é o vértice inferior-esquerdo, que é o que
satisfaz o requisito 12 do OGC Core e o que torna o EAERS requisito 27
parcialmente atendido. O centroide é conveniência.

*Rejeitado:* expor apenas o centroide. Quebraria a conformidade que o artigo
reivindica.

### D-0.5 — `nominal_cell_area` separado de `effective_cell_area`

A primeira recebe resolução, a segunda recebe célula. Assinaturas diferentes
forçam o chamador a perceber que são coisas diferentes. Células trapezoidais têm
área efetiva ≠ nominal, e num sistema cadastral retornar a nominal para uma
célula clipada é erro grave.

*Rejeitado:* uma função só com flag.

### D-0.6 — Ambos os blobs binários entram no pacote

TreeBlob **e** GeometryBlob. O TreeBlob generaliza para qualquer DGGS com índice
prefixado; o GeometryBlob é ancorado ao ITACaRT e carrega densificação e
canonização, que são geometria pura, não lógica de negócio.

*Rejeitado:* deixar o GeometryBlob fora por ser "específico do app". A identidade
jurídica on-chain é uma *aplicação* do formato, não o formato.

### D-0.7 — Aliases compatíveis com H3 v4

`latlng_to_cell`, `cell_to_latlng`, `cell_to_parent`, `cell_to_children`, com
ordem de argumentos `(lat, lng)` do H3. Reduz a barreira de adoção para quem já
usa H3.

### D-0.8 — Python ≥ 3.10

O `itacart_core` já é 3.11+. A sintaxe `X | Y` em anotações e `Literal` sem
`typing_extensions` pedem ≥ 3.10.

*Rejeitado:* manter 3.8 do placeholder.

## 5. Bugs conhecidos do protótipo

Herdados do notebook. Cada um vira teste de regressão na fase que tocar o módulo.

| ID | Descrição | Fase |
|---|---|---|
| `B-0.1` | Classes de resolução par/ímpar **invertidas** no parser: pares devem ser 1-para-4 (quaternário), ímpares 1-para-25 (quinário) | F2 |
| `B-0.2` | `get_sinusoidal_main_coordinates` recursivo — O(n²) na profundidade | F3 |
| `B-0.3` | Espelhamento de quadrante incompleto: `abs()` + reaplicação de sinal não generaliza para SW e NW | F3 |
| `B-0.4` | Sem `index_to_latlon`: existe coordenada→índice, falta o caminho de volta | F3 |
| `B-0.5` | Sem implementação de vizinhança | F6 |

O `constants.py` já carrega um `WARNING` na docstring de `QUINARY_CODES` para
que `B-0.1` não retorne.

## 6. Alvo de conformidade OGC

Frames 4 e 5 do artigo definem o que o pacote deve provar em CI.

**DGGS Core (Topic 21 / ISO 19170-1):** conformidade total, requisitos 6 a 19.

**EAERS:** conformidade parcial, e as divergências são deliberadas:

- Req 22–25 **não atendidos** — tesselação direta na superfície, sem interface
  poliédrica. É escolha de projeto, não deficiência.
- Req 21, 28–29 **parciais** — trapézios do antemeridiano não têm área igual
- Req 27 **parcial** — posição representativa é vértice, não centroide

A suite em `tests/conformance/` pina essas divergências como intencionais, para
que não virem regressão silenciosa.

## 7. Toolchain

- **Runtime:** `shapely` ≥ 2.0 (única dependência obrigatória)
- **Extras:** `geo` (GeoPandas), `parallel` (joblib), `dev`, `docs`
- **Qualidade:** `black`, `isort`, `flake8`, `mypy --strict`, `pytest` + `pytest-cov`
- **CI:** GitHub Actions — lint, matriz 3 SOs × 4 versões Python, job de
  conformidade separado, build com `twine check`
- **Docs:** Sphinx + `furo` + `napoleon` + `autodoc-typehints`, deploy no GitHub
  Pages com `-W` (warnings como erro)
- **Publicação:** PyPI via trusted publishing (OIDC, sem token)

## 8. Restrições de ambiente

Herdadas do `itacart-app`, já documentadas e com workaround conhecido.

**Windows — SSL:** o cert store pode ter certificado corrompido
(`ASN1: NOT_ENOUGH_DATA`), quebrando import de `aiohttp` (puxado por
`contextily`). O monkey-patch já está em `tests/conftest.py`. Correção
permanente, como administrador:

```powershell
certutil -generateSSTFromWU root.sst
certutil -addstore -f root root.sst
```

**Conda — numpy 2.4.0:** PEP 489 rejeita reimport. Se aparecer
`cannot load module more than once per process`, usar `numpy<2.4`.

## 9. Estilo de comunicação

- Conciso e direto, com headers
- Em bifurcação, apresentar alternativas com tradeoffs explícitos — a decisão é
  do Israel
- Não oscilar entre versões sem orientação (evitar revisões silenciosas)
- Quando ele diz "não é assim", incorporar sem defensividade
- Diagramas Mermaid **não renderizam** no fluxo dele; usar descrição textual
- Maximizar gráficos e imagens nos entregáveis (orientação do orientador no IEAv)

## 10. Contexto acadêmico

Israel é mestrando em Ciências e Tecnologias Espaciais no ITA e oficial da FAB.
Este pacote é a **implementação de referência** citada como trabalho futuro na
seção 5 do artigo, e insumo para:

- Dissertação de mestrado em curso
- **Grid Lab** — workbench de e-learning para DGGS, submissão CATCON 9 no
  XXV Congresso ISPRS (Toronto, julho de 2026)
- Benchmarking quantitativo contra H3, S2 e rHEALPix

O pacote existir e ser instalável é pré-requisito para a demonstração no
congresso.
