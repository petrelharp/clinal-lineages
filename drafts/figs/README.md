# Figures for the manuscript:

Here is the description of some scripts that produce figures in the manuscript (and/or ones similar to those).
These are often modifications of scripts in the [sims/](../../sims/) directory. These scripts will require editing to specify parameters, or an input file.

- `make_length_along_chr.R` provides examples of how to turn these into plots and heatmaps as they appear in the manuscript (e.g. Figure 5). 
    This is basically identical to `blockLengthPlotsS0.001.R`.
    This particular script will produce a bunch of figures for a simulation we ran under s=0.01, T=1000 and sigma=1, namely:

    * `blocksAlongChromAncBConditioning.pdf` (e.g. Figure 5)
    * `blocksAlongChromHeatmapAncBConditioning.pdf` (e.g. Figure S10)
    * `blocksAlongChromNoConditioning.pdf` (e.g. Figure S11)
    * `blocksAlongChromHeatmapNoConditioning.pdf` (e.g. Figure S12)
    * `adjacentBlocksAlongChromAncBConditioning.pdf` (e.g. Figure S13)
    * `ratioAdjacentBlocksAlongChromAncBConditioning.pdf` (e.g. Figure S14)
    * `ratioAdjacentBlocksAlongChromHeatmapAncBConditioning.pdf` (e.g. Figure S15)

- `Plot_frequency_clines.R` will plot ancestry frequency clines for multiple loci across the genome for both (Corresponds to Fig. 2 in the manuscript). 

- `number_of_ancestors.R` will plot the average family size of migrants. You will need to run `number_of_inds.R` first. 
    A lot of variables are hard-coded (sorry!), so these will need to be modified to suit your needs.

- `theory-hap-comparison.R` compares some haplotype-based stats from theory to simulation, producing the figure `cond_freqs_comparison_tau_<T>.pdf` (e.g. Figure S16)

- `theory-figs.R` produces 2D plots of expected allele frequency against space and time for a range of distances from the selected sites based on numerical solutions of the theory (e.g. Figure 2)

- `theory-hap-lens.R` produces some plots that aren't in the manuscript for haplotype length predictions based on theory. 
