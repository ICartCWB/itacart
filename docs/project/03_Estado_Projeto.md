# 03 — Estado do Projeto

> **Documento vivo.** Atualizado pelo chat-ponte ao fechar cada fase. É o que
> permite substituir o chat-ponte sem perder informação — se este arquivo estiver
> desatualizado, a substituição perde contexto.

**Última atualização:** F0 aberta · agosto de 2026
**Chat-ponte atual:** #1 (primeiro)

---

## Estado atual

### Resumo executivo

Esqueleto do pacote completo e versionado. Nenhuma fase de implementação
concluída. Próximo passo: F0 (bootstrap).

| Métrica | Valor |
|---|---|
| Fases concluídas | 0 de 11 |
| Funções públicas definidas | 105 |
| Funções implementadas | 0 |
| Testes | 3 passed, 15 xfailed |
| Cobertura | — |
| Versão | 0.1.0a4 (não publicada) |

### Progresso por fase

| Fase | Escopo | Status | Testes | Cobertura |
|---|---|---|---|---|
| F0 | Bootstrap, `constants`, `exceptions` | ⬜ | — | — |
| F1 | `geodesy` | ⬜ | — | — |
| F2 | `index` | ⬜ | — | — |
| F3 | `resolutions` + `cells` | ⬜ | — | — |
| F4 | `boundary` | ⬜ | — | — |
| F5 | `hierarchy` | ⬜ | — | — |
| F6 | `topology` | ⬜ | — | — |
| F7 | `geometry` | ⬜ | — | — |
| F8 | `serialization` | ⬜ | — | — |
| F9 | `interop` + `engine` + conformidade | ⬜ | — | — |
| F10 | Release v1.0.0 | ⬜ | — | — |

### O que existe no repositório

```
itacart/
├── pyproject.toml              # completo: deps, extras, black/isort/mypy/pytest
├── README.md                   # completo
├── CHANGELOG.md  CONTRIBUTING.md  .gitignore
├── .github/workflows/          # ci.yml, docs.yml, publish.yml
├── docs/
│   ├── project/                # 00 a 04 (este sistema)
│   ├── source/                 # Sphinx: conf.py, index.rst, api/, concepts/
│   └── _static/                # vazio — figuras do artigo entram aqui
├── src/itacart/
│   ├── __init__.py             # 109 nomes exportados + aliases H3
│   ├── constants.py            # STUB — F0
│   ├── exceptions.py           # STUB — F0
│   ├── geodesy.py              # STUB — F1
│   ├── index.py                # STUB — F2
│   ├── resolutions.py          # STUB — F3
│   ├── cells.py                # STUB — F3
│   ├── boundary.py             # STUB — F4
│   ├── hierarchy.py            # STUB — F5
│   ├── topology.py             # STUB — F6
│   ├── geometry.py             # STUB — F7
│   ├── serialization/          # STUB — F8
│   ├── interop.py              # STUB — F9
│   ├── engine.py               # STUB — F9
│   └── py.typed
└── tests/
    ├── conftest.py             # fixtures geográficos + workaround SSL Windows
    ├── unit/                   # um arquivo por módulo, todos vazios
    └── conformance/            # 14 testes em xfail, um por requisito OGC
```

Todos os stubs levantam `NotImplementedError`, com type hints e docstrings
completas. O contrato de API está fechado; falta o corpo.

### Bloqueios ativos

Nenhum.

---

## Registro de decisões

Formato: código, título, decisão, justificativa, alternativa rejeitada.

### Decisões de fundação (pré-F0)

Tomadas durante o desenho do esqueleto. Detalhamento completo em
`01_Briefing_Mestre.md` seção 4.

| Código | Título | Resumo |
|---|---|---|
| `D-0.1` | Escopo completo na v1 | Fronteiras e vizinhança não são adiáveis: sem elas o Req 8–10 do OGC não fecha |
| `D-0.2` | Semântica vetorizada | Índice composicional em toda a API, retorno alinhado posicionalmente. Rejeitados tipos `Cell`/`Region` separados |
| `D-0.3` | Sem `pyproj` em runtime | Eq. (1) e (2) implementadas diretamente; `pyproj` só em teste |
| `D-0.4` | `cell_to_anchor` ≠ `cell_to_centroid` | Anchor satisfaz Req 12 e explica o Req 27 parcial |
| `D-0.5` | `nominal_cell_area` ≠ `effective_cell_area` | Assinaturas diferentes forçam a distinção; trapézio tem área ≠ nominal |
| `D-0.6` | Ambos os blobs no pacote | TreeBlob generaliza, GeometryBlob é geometria ITACaRT pura |
| `D-0.7` | Aliases H3 v4 | Reduz barreira de adoção |
| `D-0.8` | Python ≥ 3.10 | `itacart_core` já é 3.11+; sintaxe de anotações pede 3.10 |

### Decisões por fase

*Vazio — nenhuma fase concluída.*

<!--
Formato ao preencher:

### F<N> — <nome da fase>

#### D-<N>.<M> — <título>

**Decisão.** O que foi decidido.
**Justificativa.** Por quê.
**Alternativa rejeitada.** O que não foi feito e por quê — evita re-litígio.
**Gatilho de reabertura.** Que condição justificaria revisitar (quando aplicável).
-->

---

## Pendências abertas

| Código | Descrição | Bloqueia | Origem |
|---|---|---|---|
| `P-0.1` | Configurar GitHub Pages no repositório (Settings → Pages → GitHub Actions) | Deploy da documentação | F0 |
| `P-0.2` | Configurar trusted publishing no PyPI para `ICartCWB/itacart` | Publicação em F10 | F0 |
| `P-0.3` | Obter figuras do artigo em alta resolução para `docs/_static/` | Páginas de conceito em F10 | F0 |
| `P-0.4` | Verificar se `P-F3.1` do `itacart-app` (bug em `geometry_to_tree`) foi resolvido lá | F8 | Handoff `itacart-app` |
| `P-0.5` | Confirmar se o patch `P-2.2.9` (§17 envelope) foi aplicado na spec do GeometryBlob | F8 | Handoff `itacart-app` |

`P-0.1` a `P-0.3` são tarefas operacionais do Israel, fora do fluxo de chats.

---

## Bugs conhecidos

Herdados do protótipo Colab. Cada um vira teste de regressão permanente na fase
que tocar o módulo.

| Código | Descrição | Fase | Status |
|---|---|---|---|
| `B-0.1` | Classes par/ímpar invertidas no parser: par é 1-para-4, ímpar é 1-para-25 | F2 | ⬜ |
| `B-0.2` | `get_sinusoidal_main_coordinates` recursivo, O(n²) na profundidade | F3 | ⬜ |
| `B-0.3` | Espelhamento de quadrante incompleto — falha em SW e NW | F3 | ⬜ |
| `B-0.4` | Ausência de `index_to_latlon` (caminho inverso) | F3 | ⬜ |
| `B-0.5` | Vizinhança não implementada | F6 | ⬜ |

### Corrigidos antes de F0

| Código | Descrição | Correção |
|---|---|---|
| `B-0.6` | `build-system` exigia `setuptools>=61`, incompatível com `license = "MIT"` em forma SPDX (PEP 639) | Elevado para `>=77`; `python -m build` e `twine check` verificados |
| `B-0.7` | `LICENSE` ausente no esqueleto, quebrando `license-files` | Arquivo MIT criado |

---

## Métricas de referência

Alvos herdados do `itacart-app`, onde já foram atingidos.

| Item | Alvo | Fonte |
|---|---|---|
| Cobertura por módulo | ≥ 85% | Disciplina do projeto |
| Cobertura `itacart_core` original | ~99% | `itacart-app` F1 |
| Concordância Vincenty vs `pyproj` | < 1 mm em 5 500 km | `itacart-app` F1+ (medido: 12 µm) |
| TreeBlob, folha resolução 13 | 10 bytes exatos | `itacart-app` F2.1 |
| TreeBlob, densidade assintótica | ~5,13 b/v em 1 000 vértices | `itacart-app` F2.1 |
| GeometryBlob, densidade | 27–38 b/v em datasets reais | `itacart-app` F3 |
| `geo_to_cell` resolução 13 | < 1 ms por ponto | Novo (F3) |

---

## Histórico de fases

*Vazio — nenhuma fase concluída.*

<!--
Formato ao preencher:

### F<N> — <nome> · <data> · ✅

**Entregue.** Lista dos módulos e artefatos.
**Testes.** N testes, X% de cobertura.
**Decisões.** Códigos D-<N>.<M> registrados.
**Pendências abertas.** Códigos P-<N>.<M>.
**Bugs corrigidos.** Códigos B-<N>.<M> com teste de regressão.
**Desvios do briefing.** O que saiu diferente do planejado e por quê.
-->

---

## Histórico de chats-ponte

| # | Período | Fases cobertas | Motivo da substituição |
|---|---|---|---|
| 1 | ago/2026 — | F0 → | — (ativo) |
