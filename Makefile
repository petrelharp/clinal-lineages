THIS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
MATHJAX = https://cdn.mathjax.org/mathjax/latest/MathJax.js
PANDOC_OPTS =  --to html --from markdown-implicit_figures --self-contained --standalone --section-divs --template /usr/local/lib/R/site-library/rmarkdown/rmd/h/default.html --variable 'theme:bootstrap' --include-in-header $(THIS_DIR)/resources/header-scripts.html --mathjax --variable 'mathjax-url:$(MATHJAX)?config=TeX-AMS-MML_HTMLorMML' --no-highlight --variable highlightjs=/usr/local/lib/R/site-library/rmarkdown/rmd/h/highlight 

%.md : %.Rmd
	cd $$(dirname $<); Rscript -e 'knitr::knit(basename("$<"),output=basename("$@"))'

%.html : %.md
	pandoc $< $(PANDOC_OPTS) --output $@

