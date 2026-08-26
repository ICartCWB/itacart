# Minimal makefile for the itacart documentation.
#
# The build imports the package from ../src, so no installation is needed.
# `make html` is the usual target; `make strict` is what CI runs.

SPHINXBUILD ?= sphinx-build
SOURCEDIR    = source
BUILDDIR     = _build

.PHONY: help html strict live coverage linkcheck clean

help:
	@echo "html      build the HTML documentation into $(BUILDDIR)/html"
	@echo "strict    build with warnings as errors, as CI does"
	@echo "live      rebuild and reload on save (needs sphinx-autobuild)"
	@echo "coverage  list public names that carry no docstring"
	@echo "linkcheck verify every external link"
	@echo "clean     remove the build output and the generated pages"

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"
	@echo
	@echo "Open $(BUILDDIR)/html/index.html"

strict:
	$(SPHINXBUILD) -b html -W --keep-going "$(SOURCEDIR)" "$(BUILDDIR)/html"

live:
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" --watch ../src

coverage:
	$(SPHINXBUILD) -b coverage "$(SOURCEDIR)" "$(BUILDDIR)/coverage"
	@cat "$(BUILDDIR)/coverage/python.txt"

linkcheck:
	$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck"

clean:
	rm -rf "$(BUILDDIR)"
	rm -rf "$(SOURCEDIR)/api/generated"
	rm -rf "$(SOURCEDIR)/_generated"
