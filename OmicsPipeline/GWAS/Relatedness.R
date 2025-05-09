library("optparse")
option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains"))
args <- parse_args(OptionParser(option_list=option_list))


pdf(paste0(args$dir,"/","relatedness.pdf"))
relatedness = read.table(paste0(args$dir,"/","pihat_min0.2.genome"), header=T)
par(pch=16, cex=1)
with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))
legend(1,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))

pdf(paste0(args$dir,"/","zoom_relatedness.pdf"))
relatedness_zoom = read.table(paste0(args$dir,"/","zoom_pihat.genome"), header=T)
par(pch=16, cex=1)
with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n"))
with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))
legend(0.02,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))

pdf(paste0(args$dir,"/","hist_relatedness.pdf"))
relatedness = read.table(paste0(args$dir,"/","pihat_min0.2.genome"), header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")  
dev.off() 

