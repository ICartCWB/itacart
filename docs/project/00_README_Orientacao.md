# 00 — Orientação do Sistema Multi-Chat

> **Leia este arquivo primeiro.** Ele explica como o desenvolvimento do pacote
> `itacart` está organizado e o que cada tipo de conversa pode ou não fazer.

---

## Por que existe este sistema

O pacote tem ~105 funções públicas distribuídas em 11 módulos com dependências
encadeadas. Implementar tudo em uma conversa só esgota o contexto antes da
metade. A divisão em fases resolve isso: cada fase é uma unidade fechada de
trabalho, executada em uma conversa própria, que termina produzindo um handoff.

O custo desse arranjo é a perda de continuidade. Os quatro documentos numerados
compensam isso — eles são a memória persistente do projeto.

---

## Os dois tipos de conversa

### Chat-ponte (orquestrador)

Coordena o projeto. **É o único autorizado a modificar os documentos 00 a 04.**

Responsabilidades:

- Abrir fases, fornecendo ao chat-filho o briefing correto
- Receber handoffs e incorporá-los em `03_Estado_Projeto.md`
- Registrar decisões arquiteturais (códigos `D-X.Y`) e pendências (`P-X.Y`)
- Decidir sequenciamento e reagir a bloqueios
- Quando o contexto encher: produzir um handoff de ponte e ser substituído

Um chat-ponte **não escreve código de produção**. Se ele começar a implementar,
o contexto enche e o papel se perde.

### Chat-filho (executor)

Executa uma fase específica. Produz código, testes e documentação técnica.

Responsabilidades:

- Implementar o escopo da fase até os critérios de aceite
- Escrever testes com cobertura ≥ 85%
- Registrar decisões tomadas durante a execução
- Terminar com um bloco de handoff no formato de `04_Template_Handoff.md`

Um chat-filho **não modifica os documentos 00 a 04**. Ele reporta; o chat-ponte
consolida.

---

## Como abrir uma fase

O chat-ponte fornece ao chat-filho, em uma conversa nova:

1. `01_Briefing_Mestre.md` — contexto do projeto inteiro
2. O bloco da fase em `02_Fases_Briefing.md`
3. A seção "Estado atual" de `03_Estado_Projeto.md`
4. O artigo ITACaRT em PDF (especificação canônica)
5. Os stubs do módulo a implementar

Prompt de abertura sugerido:

```
Você é um chat-filho executando a fase F<N> do projeto itacart.
Leia o briefing mestre, o bloco da fase e o estado atual anexados.
Implemente o escopo até os critérios de aceite. Não modifique os
documentos 00 a 04 — ao final, produza o bloco de handoff no formato
do template.
```

## Como fechar uma fase

O chat-filho produz o handoff. O chat-ponte então:

1. Verifica os critérios de aceite contra o entregue
2. Incorpora decisões e pendências em `03_Estado_Projeto.md`
3. Marca a fase como concluída em `02_Fases_Briefing.md`
4. Commita os documentos atualizados no repositório

---

## Substituição do chat-ponte

Quando o contexto do chat-ponte encher, ele produz um **handoff de ponte** — que
não é o mesmo que o handoff de fase. O handoff de ponte é apenas:

> "O estado está em `03_Estado_Projeto.md`, atualizado até a fase F<N>.
> As pendências abertas são <lista>. A próxima fase é F<N+1>."

Isso é possível porque `03_Estado_Projeto.md` é mantido completo e auto-suficiente
a cada fase fechada. **Se esse documento não estiver em dia, a substituição do
chat-ponte perde informação.** Mantê-lo atualizado é a obrigação mais importante
do papel.

---

## Disciplina de código

Herdada do projeto `itacart-app`, que já validou essas convenções:

- Docstrings e identificadores em inglês; comentários em português técnico
  quando ajudam contexto
- Type hints obrigatórios em toda assinatura pública
- Cobertura alvo ≥ 85% por módulo
- Bug corrigido ganha teste de regressão permanente
- Decisão arquitetural ganha código `D-X.Y` com justificativa **e alternativas
  rejeitadas** — as rejeições evitam re-litígio em fases futuras
- Pendência ganha código `P-X.Y`
- Toda fase termina com `black`, `isort`, `flake8`, `mypy` e `pytest` limpos

---

## Convenção de códigos

| Prefixo | Significado | Exemplo |
|---|---|---|
| `F<N>` | Fase | `F4` — Comportamento de fronteira |
| `D-<N>.<M>` | Decisão arquitetural | `D-4.1` — trapézio calculado por clipping |
| `P-<N>.<M>` | Pendência aberta | `P-2.1` — validar contra dataset real |
| `B-<N>.<M>` | Bug com teste de regressão | `B-3.2` — sinal invertido no quadrante SW |

O número da fase vem primeiro, o sequencial dentro dela depois.

---

## Fonte canônica de verdade

Quando houver divergência entre documentos, esta é a ordem de precedência:

1. **O artigo publicado** (DOI 10.14393/rbcv77n0a-79281) — especificação do DGGS
2. `01_Briefing_Mestre.md` — decisões de projeto do pacote
3. `03_Estado_Projeto.md` — estado corrente
4. Docstrings nos stubs — contrato de API
5. Handoffs individuais — histórico

Se um chat-filho encontrar contradição entre o artigo e um stub, o artigo vence
e a divergência vira um item de handoff.

---

## Arquivos deste sistema

| Arquivo | Conteúdo | Quem modifica |
|---|---|---|
| `00_README_Orientacao.md` | Este arquivo | Chat-ponte |
| `01_Briefing_Mestre.md` | Contexto e decisões do projeto | Chat-ponte |
| `02_Fases_Briefing.md` | Briefing de cada fase | Chat-ponte |
| `03_Estado_Projeto.md` | Estado corrente, decisões, pendências | Chat-ponte |
| `04_Template_Handoff.md` | Formato do handoff de fase | Chat-ponte |

Todos versionados em `docs/project/` no repositório.
