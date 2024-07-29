# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
PACKAGE_DIR   = ../src

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Custom targets for generating API docs
.PHONY: apidoc
apidoc: apidoc-public apidoc-private

apidoc-public:
	sphinx-apidoc --templatedir=source/_templates/public -f -o "$(SOURCEDIR)/public_api" "$(PACKAGE_DIR)"

apidoc-private:
	sphinx-apidoc --private --templatedir=source/_templates/private -f -o "$(SOURCEDIR)/private_api/" "$(PACKAGE_DIR)"

# Python script to customize .rst files
.PHONY: customize-rst
customize-rst:
	python rename.py

# Modify the html target to depend on the new apidoc and customize-rst targets
html: apidoc customize-rst
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	@echo "Cleaning up build directories..."
	rm -r "$(BUILDDIR)"
	rm -r "$(SOURCEDIR)/auto_examples"
	rm "$(SOURCEDIR)/_static/"*
	@echo "Clean-up done."


# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
