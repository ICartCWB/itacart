# Contribuindo

## Ambiente

```bash
git clone https://github.com/ICartCWB/itacart
cd itacart
pip install -e ".[dev,docs,geo]"
```

## Antes de abrir PR

```bash
black src tests
isort src tests
flake8 src tests --max-line-length=88 --extend-ignore=E203,W503
mypy src
pytest
```

## Disciplina de codigo

- Docstrings e identificadores em ingles; comentarios em portugues tecnico quando ajudam contexto.
- Type hints obrigatorios em toda assinatura publica.
- Cobertura alvo >= 85% por modulo.
- Bug corrigido ganha teste de regressao permanente.
- Decisao arquitetural relevante entra em `docs/project/03_Estado_Projeto.md` com justificativa e alternativas rejeitadas — inclusive as rejeicoes, para evitar re-litigio.

## Ordem sugerida de implementacao

O grafo de dependencias entre modulos define a ordem:

1. `constants`, `exceptions` — sem dependencias
2. `geodesy` — projecao e Vincenty
3. `index` — parser, compose, normalize
4. `resolutions` — tabela e escalas
5. `cells` — quantizacao e geometria inversa
6. `boundary` — fronteiras (depende de cells)
7. `hierarchy` — navegacao
8. `topology` — vizinhanca (depende de boundary para deflexao)
9. `geometry` — polyfill, densificacao, canonizacao
10. `serialization` — TreeBlob e GeometryBlob
11. `interop`, `engine`

## Sistema multi-chat

O desenvolvimento e coordenado por fases documentadas em `docs/project/`.
Ver `docs/project/00_README_Orientacao.md`.
