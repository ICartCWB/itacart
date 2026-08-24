# 04 — Template de Handoff

> **Para o chat-filho.** Ao terminar a fase, produza um bloco no formato abaixo.
> O chat-ponte consome esse bloco para atualizar `03_Estado_Projeto.md`.
>
> Preencha todas as seções. Uma seção vazia deve dizer "nenhum", não ser omitida —
> a ausência é ambígua, o "nenhum" é informação.

---

## Formato

```markdown
# Handoff — F<N>: <nome da fase>

**Data:** <mês/ano>
**Status:** ✅ concluída | ⚠️ parcial | ❌ bloqueada

## 1. Entregue

Módulos implementados e artefatos produzidos. Caminhos relativos à raiz do repo.

- `src/itacart/<modulo>.py` — <N> funções, <M> linhas
- `tests/unit/test_<modulo>.py` — <N> testes
- <outros artefatos: specs, fixtures, notebooks>

## 2. Critérios de aceite

Um item por critério do briefing da fase, com verificação explícita.

| # | Critério | Status | Evidência |
|---|---|---|---|
| 1 | <texto resumido do critério> | ✅ | `test_nome_do_teste` |
| 2 | <...> | ⚠️ | <o que faltou e por quê> |

## 3. Testes e cobertura

- **Total:** <N> testes, <N> passed, <N> xfailed, <N> skipped
- **Cobertura:** <X>% no módulo, <Y>% no pacote
- **Tempo de execução:** <T>s
- **Lint:** `black` ✅ · `isort` ✅ · `flake8` ✅ · `mypy` ✅

Se algum estiver vermelho, explicar.

## 4. Decisões tomadas

Uma entrada por decisão que alguém poderia questionar depois.

### D-<N>.<M> — <título>

**Decisão.** O que foi decidido.

**Justificativa.** Por que essa e não outra.

**Alternativa rejeitada.** O que foi considerado e descartado, com o motivo.
Registrar a rejeição evita que a fase seguinte re-litigue a mesma questão.

**Gatilho de reabertura.** Que condição futura justificaria revisitar.
Omitir quando não houver.

*Se nenhuma decisão relevante foi tomada, escrever "Nenhuma."*

## 5. Bugs corrigidos

| Código | Descrição | Teste de regressão |
|---|---|---|
| `B-<N>.<M>` | <o que estava errado> | `test_<nome>` |

Todo bug corrigido ganha teste permanente. Sem teste, o bug não conta como
corrigido.

## 6. Pendências abertas

| Código | Descrição | Bloqueia | Severidade |
|---|---|---|---|
| `P-<N>.<M>` | <o que ficou pendente> | <fase ou nada> | alta/média/baixa |

Distinguir pendência (trabalho conhecido não feito) de bug (comportamento
errado).

## 7. Desvios do briefing

O que saiu diferente do planejado, e por quê. Inclui escopo que cresceu,
encolheu, ou mudou de forma.

*Se a fase seguiu o briefing, escrever "Nenhum."*

## 8. Notas para a próxima fase

O que quem pegar a fase seguinte precisa saber e não está óbvio no código.
Armadilhas encontradas, suposições feitas, partes frágeis.

## 9. Estado do repositório

- **Branch:** <nome>
- **Commits:** <N> commits nesta fase
- **CI:** ✅ verde | ❌ vermelho (<motivo>)
- **Arquivos modificados fora do escopo da fase:** <lista ou "nenhum">

Modificar arquivo fora do escopo não é proibido, mas precisa ser declarado —
é como o chat-ponte detecta acoplamento não previsto.
```

---

## Notas sobre o preenchimento

**Seja específico em "Evidência".** "Testado" não é evidência; `test_req_13_unique_address`
é. O chat-ponte não tem acesso ao seu contexto e precisa poder verificar.

**Decisões rejeitadas são o item mais valioso.** A fase seguinte não sabe o que
você já considerou e descartou. Sem esse registro, ela reabre a discussão do zero
ou, pior, toma a decisão oposta sem perceber.

**Desvios não são falha.** O briefing foi escrito antes de a fase começar; que
ele erre em algum ponto é esperado. O que causa problema é o desvio não
declarado, que só aparece três fases depois.

**Se a fase estourou o contexto**, produza o handoff mesmo assim, com status
⚠️ parcial, listando exatamente o que ficou de fora. Um handoff parcial preciso
vale mais que um completo impreciso — o chat-ponte abre uma fase de continuação
a partir dele.

---

## Exemplo abreviado

```markdown
# Handoff — F2: Índice composicional

**Data:** set/2026
**Status:** ✅ concluída

## 1. Entregue

- `src/itacart/index.py` — 11 funções, 340 linhas
- `tests/unit/test_index.py` — 62 testes

## 2. Critérios de aceite

| # | Critério | Status | Evidência |
|---|---|---|---|
| 1 | Figura 7 faz round-trip sem perda | ✅ | `test_figure_7_roundtrip` |
| 5 | B-0.1 corrigido (par/ímpar) | ✅ | `test_even_rejects_quinary_codes` |

## 4. Decisões tomadas

### D-2.1 — `normalize` ordena irmãos pelo alfabeto de refinamento

**Decisão.** Ordenação pelo índice no alfabeto (`1`<`2`<`3`<`4`;
`A1`<`A2`<...<`E5`), não lexicográfica sobre a string.

**Justificativa.** Ordem lexicográfica colocaria `A10` antes de `A2` se o
alfabeto crescesse. Ordenar pela posição no alfabeto é estável a extensões.

**Alternativa rejeitada.** `sorted()` sobre a string. Funciona hoje porque o
alfabeto vai só até 5, mas é frágil.

## 8. Notas para a próxima fase

`decompose` é O(n) mas materializa a lista inteira. F7 (`polyfill`) em resolução
13 deve usar `iter_cells`, que é gerador.
```
