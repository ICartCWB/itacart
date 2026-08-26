# Building this documentation

The build imports the package from `src/`, not from `site-packages`. **No
installation is required**, and an installed copy that shadows the working tree
makes `conf.py` fail loudly rather than document the wrong code.

## Once

```console
$ pip install -e ".[docs]"
```

The `-e` is optional — it installs the package, but the build does not use the
installed copy. What the extra provides is Sphinx itself, the Furo theme,
MyST and the two Sphinx plugins.

## Build

```console
$ cd docs
$ make html
```

The result is `docs/_build/html/index.html`; open it in a browser.

On Windows, `make.bat html` from the same directory.

## Live preview

`sphinx-autobuild` rebuilds on save and reloads the browser, which is the
useful mode while writing docstrings:

```console
$ sphinx-autobuild docs/source docs/_build/html --watch src
```

The `--watch src` is the point: a docstring edited in `src/itacart/` rebuilds
the page that documents it.

## What is generated, and what is written

Almost nothing in `docs/source/` is maintained by hand.

| Path | Origin |
|---|---|
| `api/generated/` | `sphinx.ext.autosummary`, walking the package recursively |
| `_generated/figures/` | the `phase_figures` extension, scanning `docs/_static/fN/` |
| `_build/` | Sphinx output |
| everything else | written by hand: `conf.py`, `index.md`, this page, `api/index.md` |

Both generated directories are rewritten on every build and should be ignored
by Git.

Adding a module in a later phase requires **no edit here**: `autosummary`
finds it through `itacart`'s own namespace. Adding a phase's figures requires
no edit either — drop the PNGs in `docs/_static/f5/` and the gallery page
appears, with captions derived from the filenames.

## Checking for undocumented code

```console
$ cd docs
$ make coverage
$ cat _build/coverage/python.txt
```

`sphinx.ext.coverage` lists every public name reached by autodoc that carries
no docstring. On a package whose docstrings are the specification, that list
should stay empty.

## Strict builds

The CI job builds with warnings as errors, and so should a local check before
pushing:

```console
$ sphinx-build -b html -W --keep-going docs/source docs/_build/html
```

A broken cross-reference, a missing image or a malformed docstring becomes a
build failure rather than a silently ugly page.
