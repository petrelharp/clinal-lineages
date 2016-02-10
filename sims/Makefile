include ../Makefile

%_comparison-to-theory.md : %_simsums_chunks.Robj
	Rscript -e "library(templater); simchunks.file='$<'; render_template('compare-to-theory.Rmd',output='$@')"

%_comparison-to-theory.html : %_simsums_chunks.Robj
	Rscript -e "library(templater); simchunks.file='$<'; render_template('compare-to-theory.Rmd',output='$@')"
