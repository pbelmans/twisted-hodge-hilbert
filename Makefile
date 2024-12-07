SPHINXOPTS    ?=
SPHINXBUILD   ?= sage --python -msphinx
SOURCEDIR     =
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


install:
	sage -pip install --upgrade -v -e .
	rm -rf twisted_hodge_hilbert.egg-info

test:
	sage -t twisted_hilbert

coverage:
	sage --coverage twisted_hilbert

lint:
	mkinit --black twisted_hilbert/__init__.py > twisted_hilbert/__init__.py
	black twisted_hilbert
	isort --profile black twisted_hilbert
	flake8 --extend-ignore=E741 --max-line-length 89 twisted_hilbert
	ruff check --ignore=E741 twisted_hilbert

.PHONY: install test coverage lint
