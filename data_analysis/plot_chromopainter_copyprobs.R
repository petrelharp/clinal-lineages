#plot chromopainter copyprob plot.

setwd("/hdd/HZ_selection/mouse_data_HZ/Unpruned_mouse_data/unpruned_mouse_by_chrom/")
files =list.files()
probfiles = files[which(grepl("unpruned",files)==F)]

files_inds = sapply(probfiles,function(X){paste("X_",strsplit(strsplit(X,"_")[[1]][2],"ind")[[1]][1],sep="")})


indtable = read.table("/hdd/HZ_selection/mouse_data_HZ/mouse_by_chrom/mousedata_indfile",as.is=T)
inds = indtable$V1

pdf(file="/hdd/HZ_selection/mouse_data_HZ/Unpruned_mouse_data/chromopaintings.pdf",width=13,height=3)

for(i in 1:length(probfiles)){

	data=read.table(probfiles[i],as.is=T)
	plot(data$V1,data$V2,type="l",main=paste("chr",probfiles[i],"_pop",indtable$V2[which(inds==files_inds[i])],sep=""),ylab="Pr(ancA)",xlab="position (bp)")

	with(data,polygon(c(V1,V1[length(V1)],V1[1]), c(V2,0,0),col="red" ))
}

dev.off()

