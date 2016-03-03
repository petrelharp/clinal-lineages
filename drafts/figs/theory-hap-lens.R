# Get fine-scale haplotype length predictions.

source("../../sims/sim-fns.R",chdir=TRUE)
library(ReacTran)
library(Matrix)
source("../../notes/tran.1D.R",chdir=TRUE) # fix for log.A
source("../../notes/cline-fns.R",chdir=TRUE)
source("../../sims/chunks_fns.R",chdir=TRUE)

theory.param.list <- list(
          list( 
              sigma=1L,
              s=0.1,
              tau=1000L,
              density=500L,
              ndemes=50L,
              sim.files=file.path( "../../sims",c(
                      "simulation_SIGMA1_Ninds25000_ndemes50_s0.1_dir/tau1000_dir/results_runid_20160130_intervalSizes.Robj",
                      "simulation_SIGMA1_Ninds25000_ndemes50_s0.1_dir/tau1000_dir/results_runid_20160130_intervalSizes_allAncs.Robj"))
             ),
          list( 
              sigma=1L,
              s=0.01,
              tau=1000L,
              density=500L,
              ndemes=50L,
              sim.files=file.path( "../../sims",c(
                      "simulation_SIGMA1_Ninds25000_ndemes50_s0.01_dir/tau1000_dir/results_runid_20160127_intervalSizes.Robj",
                      "simulation_SIGMA1_Ninds25000_ndemes50_s0.01_dir/tau1000_dir/results_runid_20160127_intervalSizes_allAncs.Robj" ))
              )
      )


theory.params <- theory.param.list[[1]]
theory.params$ninds <- theory.params$ndemes * theory.params$density

filebase <- sprintf("haplotypes_SIGMA%d_density%d_ndemes%d_s%0.3f_tau%d",
                 theory.params$sigma, theory.params$density, theory.params$ndemes, theory.params$s, theory.params$tau )
robj.file <- file.path("cache",paste0(filebase,".Robj"))

if (!file.exists(robj.file)) {
    tt <- seq(0,theory.params$tau,length.out=21)
    xgrid <- setup.grid.1D(x.up=-theory.params$ndemes/2, x.down=theory.params$ndemes/2, N=50)
    fgrid <- extend_grid(xgrid)
    fwds.soln <- forwards_pde(s=theory.params$s, times=tt, grid=fgrid, sigma=theory.params$sigma )

    # do it on a fine scale: for tau=1000, looks to be over +/- 5 or 6 cM
    # more generally +/- 2/sqrt(tau)
    rgrid <- setup.grid.1D(x.up=-2/sqrt(theory.params$tau), x.down=2/sqrt(theory.params$tau), N=61)
    rgrid$x.mid <- round(rgrid$x.mid,digits=8)
    hap.soln <- forwards_backwards_haplotypes(s=theory.params$s, times=tt, xgrid=xgrid, rgrid=rgrid,
                               sigma=theory.params$sigma, 
                               fwds.grid=fgrid, fwds.soln=fwds.soln )

    dir.create("cache",showWarnings=FALSE)
    save(filebase, tt, xgrid, fgrid, fwds.soln, rgrid, hap.soln, file=robj.file)
} else {
    load(robj.file)
}

####
# plot stuff

# expected from theory  (see compare-to-theory.Rmd)
hap.pdf <- hap_cdf_to_pdf( hap.soln )

hap.len.mat <- sapply( rgrid$x.mid, function (rr) { 
               haplens <- mean_hap_len( hap.pdf, loc=rr )
               haplens[length(tt),1+seq_along(xgrid$x.mid)]
           } )

# if (FALSE) 
{

# from simulation
for (x in theory.params$sim.files){ load(x) }
positions=seq(0.48,0.52,0.001)
deme_ID = rep( 1:theory.params$ndemes, each=theory.params$density )
xx = (1:theory.params$ndemes)-0.5-theory.params$ndemes/2

# this is a (locus x deme) matrix giving mean lengths
#   for haplotype chunks of B ancestry
mean_per_deme_AncB = do.call(rbind, lapply(intervalSizes,function(P){tapply(1:length(deme_ID),deme_ID,function(Z){all_sites = P[Z,]; mean(all_sites[which(all_sites>0)],na.rm=TRUE)})}))

.spatial.legend <- function () {
    # do the legend for a matplot where one line correponds to a location
    these.lines <- unique(floor(seq(1,length(xx),length.out=10)))
    legend("topright",
           title="distance from center:",
           legend = xx[these.lines],
           cex=0.7,col=rainbow(length(xx))[these.lines],lty=1,xjust=1)
}


# this plots mean block length, for only ancestry B
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromAncBConditioning_nonnormalized_comparison.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2.5,1,0))
# simulation
matplot(100*(positions-0.5), 100*mean_per_deme_AncB, type='l', lty=1, col=rainbow(theory.params$ndemes),
        ylim=range(100*mean_per_deme_AncB,100*t(hap.len.mat),finite=TRUE),
        xlab="Distance from selected locus (cM)",
        ylab="Mean block length (cM)",
        main="Mean block length, ancestry B")
# theory
matlines(rgrid$x.mid*100, 100*t(hap.len.mat),
        lty=2, col=rainbow(length(xgrid$x.mid)) )

.spatial.legend()
dev.off()

# this plots mean block length, for only ancestry B
# and divided by the mean length for the deme across loci
pdf(height=4, width=6.25,file=paste(outstring,"blocksAlongChromAncBConditioning_comparison.pdf",sep="_"))
par(mar=c(3.5,3.5,1.5,0.5),mgp=c(2.5,1,0))
# simulation
matplot(100*(positions-0.5), sweep(mean_per_deme_AncB,2,colMeans(mean_per_deme_AncB,na.rm=TRUE),"/"), 
        type='l', lty=1, col=rainbow(theory.params$ndemes),
        ylim=c(0.5,3),
        xlab="Distance from selected locus (cM)",
        ylab="Mean block length (cM)",
        main="Mean block length, ancestry B")
# theory
matlines(rgrid$x.mid*100, t(sweep(hap.len.mat,1,rowMeans(hap.len.mat),"/")),
        lty=2, col=rainbow(length(xgrid$x.mid)) )

.spatial.legend()
dev.off()



pdf(height=4, width=6.25,file=paste("blocksAlongChromAncBConditioning_comparison.pdf",sep="_"))

par(mar=c(3.5,3.5,0.5,0.5))
plot( (positions-0.5)*100, mean_per_deme_AncB[,1]/mean(mean_per_deme_AncB[,1],na.rm=T), 
     ylim=c(0.5,3),type="l",lty=1,col="white",main="",cex.main=0.8,xlab="",ylab="")
for(i in 1:theory.params$ndemes){
    lines((positions-0.5)*100,mean_per_deme_AncB[,i]/mean(mean_per_deme_AncB[,i],na.rm=T),
          col=rainbow(theory.params$ndemes)[i])
}
mtext("Distance from selected locus (cM)",side=1,line=2.5)
mtext("Normalized mean block length",side=2,line=2.5)

# theory
matlines(rgrid$x.mid*100, t(sweep(hap.len.mat,1,rowMeans(hap.len.mat),"/")), 
        lty=1, col=rainbow(length(xgrid$x.mid)) )

.spatial.legend()

dev.off()

}

# if (FALSE) 
{
    pdf(file=paste0(filebase,"_haplotype_lens.pdf"),width=12,height=4,pointsize=10)
    layout(t(1:3))

    for (k in 1:length(tt)) {
        hap.len.mat <- sapply( rgrid$x.mid, function (rr) { 
                       haplens <- mean_hap_len( hap.pdf, loc=rr )
                       haplens[k,1+seq_along(xgrid$x.mid)]
                   } )

        # mean length against space
        matplot(rev(xgrid$x.mid), hap.len.mat, type='l', lty=1,
                main=sprintf("mean B haplotype length, t=%d",tt[k]), 
                col=rainbow(length(rgrid$x.mid)),
                xlab="spatial position", ylab="mean B haplotype length (incl. zeros)")
        k.legend <- floor(seq(1,ncol(hap.len.mat),length.out=8))
        legend("topleft", lty=1, col=rainbow(length(rgrid$x.mid))[k.legend],
                legend=sprintf("r=%0.0f cM",100*rgrid$x.mid[k.legend]) )

        # mean length along the genome
        matplot(rgrid$x.mid*100, t(hap.len.mat), type='l', lty=1,
                main=sprintf("mean B haplotype length, t=%d",tt[k]), 
                col=rainbow(length(xgrid$x.mid)),
                xlab="position along genome (cM)", 
                ylab="mean B haplotype length (incl. zeros)")
        k.legend <- floor(seq(1,nrow(hap.len.mat),length.out=8))
        legend("topleft", lty=1, col=rainbow(length(xgrid$x.mid))[k.legend],
                legend=sprintf("x=%0.2f",xgrid$x.mid[k.legend]) )

        # mean length along the genome, normalized
        matplot(rgrid$x.mid*100, t(sweep(hap.len.mat,1,rowMeans(hap.len.mat),"/")), 
                type='l', lty=1,
                main=sprintf("mean B haplotype length, normalized by location, t=%d",tt[k]), 
                col=rainbow(length(xgrid$x.mid)),
                xlab="position along genome (cM)", 
                ylab="relative mean B haplotype length (incl. zeros)")
        k.legend <- floor(seq(1,nrow(hap.len.mat),length.out=8))
        legend("topleft", lty=1, col=rainbow(length(xgrid$x.mid))[k.legend],
                legend=sprintf("x=%0.2f",xgrid$x.mid[k.legend]) )
    }
    dev.off()
}
