
dirs="tau50_dir tau100_dir tau150_dir tau250_dir tau750_dir tau1000_dir"

figs="adjacentBlocksAlongChromAncBConditioning.pdf adjacentBlocksAlongChromNoConditioning.pdf blocksAlongChromAncBConditioning.pdf blocksAlongChromHeatmapAncBConditioning.pdf blocksAlongChromHeatmap.pdf blocksAlongChromNoConditioning.pdf ratioAdjacentBlocksAlongChromAncBConditioning.pdf ratioAdjacentBlocksAlongChromHeatmapAncBConditioning.pdf ratioAdjacentBlocksAlongChromHeatmapNoConditioning.pdf ratioAdjacentBlocksAlongChromNoConditioning.pdf density.pdf ecdf.pdf mean.pdf freqplot.pdf LD.pdf adjacentBlocksAlongChromAncBConditioning_nonnormalized.pdf adjacentBlocksAlongChromNoConditioning_nonnormalized.pdf blocksAlongChromAncBConditioning_nonnormalized.pdf blocksAlongChromNoConditioning_nonnormalized.pdf ratioAdjacentBlocksAlongChromAncBConditioning_nonnormalized.pdf ratioAdjacentBlocksAlongChromNoConditioning_nonnormalized.pdf"

for figname in $figs
do
    pdfs=$((for dirz in $dirs; do find $dirz -name "*$figname"; done)|tr '\n' ' ')
    pdfjoin --outfile $figname $pdfs
done
