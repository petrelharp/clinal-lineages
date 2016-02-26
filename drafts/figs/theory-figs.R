library(ReacTran)
library(Matrix)
library(RColorBrewer)
source("../../notes/tran.1D.R",chdir=TRUE) # fix for log.A
source("../../notes/cline-fns.R",chdir=TRUE)


###
# styles

prop.pal <- colorRampPalette(c(rev(brewer.pal(n=6,name="Blues")),brewer.pal(n=6,name="Reds")))
# prop.pal <- colorRampPalette(rev(brewer.pal(n=6,name="Blues")))

#####
## from pde-clines.Rmd

## Time for the cline to get set up
sigma <- 1
xgrid <- setup.grid.1D(x.up=-floor(30*sigma), x.down=floor(30*sigma), N = 500)

r <- 0.1
s <- 0.02
zone.age <- 150  # this is tau
tt <- seq(0,zone.age,length.out=500)

# selected frequency
fwds.soln <- forwards_pde(s=s, times=tt, grid=extend_grid(xgrid), sigma=sigma)

# linked frequencies
rr <- c(0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5)
linked.solns <- lapply( rr, function (r) {
            forwards_backwards_pde(s=s, times=tt, grid=xgrid, 
                    r=r, sigma=sigma,
                    log.p=TRUE, fwds.soln=fwds.soln)
        } )
linked.clines <- lapply( linked.solns, cline_from_soln, grid=xgrid )

.imgplot <- function (x,which=1,
        xlab="generation", ylab="space", ...) {
    xmat <- x[,1+(which-1)*attr(x,"dimens")+1:attr(x,"dimens")]
    graphics::image( xmat, zlim=c(0,1), 
            xaxt='n', yaxt='n',
            ylim=c(0.3,0.7),
            ylab=ylab, xlab=xlab,
            col=prop.pal(100), 
           ... )
    graphics::contour( xmat, add=TRUE )
    add_grid_axis(attr(x,"grid"))
    add_grid_axis(side=1,x=x[,1])
}

# from deSolve:::drawlegend
.drawlegend <- function ( zlim=c(0,1), colfn=prop.pal, horizontal=FALSE ) {
    iy <- 1
    minz <- zlim[1]
    maxz <- zlim[2]
    binwidth <- (maxz - minz)/64
    ix <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iz <- matrix(ix, nrow = length(ix), ncol = length(iy))
    if (horizontal) {
        opar <- par(mar=c(par("mar"),3,0.5)[c(5,2,6,4)])
        on.exit(par(opar),add=TRUE)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = colfn(64))
        do.call("axis", list(side = 1, mgp = c(3, 1, 0), las = 1))
    } else {
        opar <- par(mar=c(par("mar"),3,0.5)[c(1,5,3,6)])
        on.exit(par(opar),add=TRUE)
        image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = colfn(64))
        do.call("axis", list(side = 2, mgp = c(3, 1, 0), las = 2))
    }
}


pdf( file="linked-frequencies.pdf", width=6.5, height= 4.2, pointsize=10 )
par(mar=c(3.5,3,2,0)+.1, mgp=c(2.3,1,0))
layout(matrix(c(1:9,9),nrow=2),widths=c(rep(1,4),0.4))
.imgplot( fwds.soln, main=expression(p(x,t)) )
.imgplot( linked.clines[[match(0.5,rr)]], main=expression(q(x,t,r==0.5)) )
for (k in c(2,4,8)) {
    .imgplot( linked.solns[[k]], which=1, 
            main=as.expression(substitute( q[A](x,t,r==thisr), list(thisr=rr[k]))) )
    .imgplot( linked.solns[[k]], which=2,
            main=as.expression(substitute( q[B](x,t,r==thisr), list(thisr=rr[k]))) )
}
.drawlegend()
dev.off()
