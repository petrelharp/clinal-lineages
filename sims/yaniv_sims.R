#setwd("~/Desktop/swarmz")
#library(multicore)
# MAKE STARTING SAMPLES

source("sim-fns.R")



qtl=list(chr1=data.frame(traits=c("underdominant"), s = c(0.1),pos=c(0.5)))


#sigma = n*sqrt(s)/2/sqrt(tau). Define the width of the range by the number of demes. sigma/sqrt(s) is the width of the cline. sqrt(tau) is rate at which cline width is increasing (should also depend on rqtl=list(data.frame(traits=c("underdominant"), s = c(0.1),pos=c(0.5)))

sims.s0.1 = lapply(c(50),sim.zone,n.ind=2000,n.deme=200,sigma=2.5)

sims.sums = lapply(sims.s0.1,spBreaks)
#save(sims.sums,file="~/Projects/HybridZones/ClineProjects/sims_sums_s0.1.Robj")




loci = seq(0.5,1,0.1)

#MAKE file.snp
    info = data.frame(cbind(1,round(loci*1000000),0,loci))
    colnames(info) = NULL
    #write.csv(info,row.names=F,file="file.snp")


#GENOTYPE ALL MY INDS
my.genos = lapply(sims.sums,function(Z){data.frame(do.call(cbind,lapply(Z$ind.ancest, geno.ind, loci)))})
   #   colnames(my.genos) = NULL

#MAKE OUR FILE [NOTE THIS TAKES A LONG TIME AND IT MIGHT BE SMARTER TO WRITE THE THING ABOVE TO FILE ONE LINE AT A TIME]
   # write.table(my.genos, file = "file.geno",row.names = FALSE,sep="")


genotypes = lapply(my.genos,function(Z){(Z[,seq(1,ncol(Z),2)]+ Z[,seq(2,ncol(Z),2)])/2})

inds.per.deme = 10
freqs = lapply(genotypes,function(X){apply(X,1,function(Z){tapply(Z,cut(1:ncol(X),breaks=seq(0,ncol(X),inds.per.deme)),mean)})})

######COMPARE TO THEORY:

xx <- seq(-5,5,length.out=500)
plot(-xx, pcline(xx,r=1), type='l', xlab="distance from cline center", ylab="prob of inherited from the left" )
rr=(loci-0.5)*SIGMA/sqrt(0.1)
#rr <- seq(0,2,length.out=11)
# pp <- lapply(c(5,10,50),function(T){sapply( rr, function (r) { pcline(xx,r=r,sigma=2.5,tau=T) } )})

SIGMA = sqrt(sum(MIG*(1:length(MIG))^2))
xx <- (1:ndemes)-0.5-ndemes/2
pp <- pcline(xx,r=rr,tau=50,sigma=SIGMA)

par(mfrow=c(1,3))

for(TIME in 1:3){
matplot(-xx,pp[[TIME]],type='l',lty=1,col=rainbow(length(rr)))
matpoints(freqs[[TIME]],x=(-100:99)+0.5,lty=2,col=rainbow(length(rr)),type="l")
}
