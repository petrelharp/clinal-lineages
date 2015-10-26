THIS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
LOCAL_MATHJAX = /home/usr-share/javascript/mathjax/MathJax.js
ifeq ($(wildcard $(LOCAL_MATHJAX)),)
	MATHJAX = https://cdn.mathjax.org/mathjax/latest/MathJax.js
else
	MATHJAX = $(LOCAL_MATHJAX)
endif
PANDOC_OPTS =  --to html --from markdown-implicit_figures --self-contained --standalone --section-divs --template $(THIS_DIR)/resources/rmarkdown-template.html --variable 'theme:bootstrap' --include-in-header $(THIS_DIR)/resources/header-scripts.html --mathjax --variable 'mathjax-url:$(MATHJAX)?config=TeX-AMS-MML_HTMLorMML' --variable 'libraries-url:$(THIS_DIR)/resources' --no-highlight --variable highlightjs=$(THIS_DIR)/resources/highlight 


.PHONY : test

%.md : %.Rmd
	cd $(dir $<) && Rscript -e 'knitr::opts_chunk$$set(fig.path=file.path("figure","$*",""),cache.path=file.path("cache","$*",""));knitr::knit(basename("$<"),output=basename("$@"))'

%.html : %.md
	cd $(dir $<) && pandoc $(notdir $<) $(PANDOC_OPTS) --output $(notdir $@)

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
