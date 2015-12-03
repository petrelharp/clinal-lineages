SHELL = /bin/bash

.PHONY : test publish

THIS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
LOCAL_MATHJAX = /usr/share/javascript/mathjax/MathJax.js
ifeq ($(wildcard $(LOCAL_MATHJAX)),)
	MATHJAX = https://cdn.mathjax.org/mathjax/latest/MathJax.js
else
	MATHJAX = $(LOCAL_MATHJAX)
endif

# keep latex macros in this file
LATEX_MACROS = macros.tex

PANDOC_HTML_OPTS =  --to html --from markdown-implicit_figures --self-contained --standalone --section-divs --template $(THIS_DIR)/resources/rmarkdown-template.html --variable 'theme:bootstrap' --include-in-header $(THIS_DIR)/resources/header-scripts.html --mathjax --variable 'mathjax-url:$(MATHJAX)?config=TeX-AMS-MML_HTMLorMML' --variable 'libraries-url:$(THIS_DIR)/resources' --no-highlight --variable highlightjs=$(THIS_DIR)/resources/highlight 
PANDOC_PDF_OPTS = 

ifeq ($(wildcard $(LATEX_MACROS)),)
	# LATEX_MACROS doesn't exist
else
	PANDOC_HTML_OPTS += -H <(echo '\['; cat $(LATEX_MACROS); echo '\]')
	PANDOC_PDF_OPTS += -H $(LATEX_MACROS)
endif



%.md : %.Rmd
	cd $(dir $<) && Rscript -e 'knitr::opts_chunk$$set(fig.path=file.path("figure","$*",""),cache.path=file.path("cache","$*",""));knitr::knit(basename("$<"),output=basename("$@"))'

%.html : %.md
	cd $(dir $<) && pandoc $(notdir $<) $(PANDOC_HTML_OPTS) --output $(notdir $@)

%.pdf : %.md
	cd $(dir $<) && pandoc $(notdir $<) $(PANDOC_PDF_OPTS) --output $(notdir $@)

# save inkscape svg files as .ink.svg and this'll do the right thing
%.svg : %.ink.svg
	inkscape $< --export-plain-svg=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

test : 
	echo "Directory of this makefile: $(THIS_DIR) ."

# update html in the gh-pages branch
#   add e.g. 'pdfs' to the next line to also make pdfs available there
#
# hope your head isn't detached
GITBRANCH := $(shell git symbolic-ref -q --short HEAD)

publish :
	@if ! git diff-index --quiet HEAD --; then echo "Commit changes first."; exit 1; fi
	-mkdir htmls
	cp $$(find . -name '*html') htmls
	git checkout gh-pages
	cp -r htmls/* .
	git add *.html
	git commit -a -m 'automatic update of html'
	git checkout $(GITBRANCH)


