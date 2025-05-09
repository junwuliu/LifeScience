library("optparse")
option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains"))
args <- parse_args(OptionParser(option_list=option_list))

hwe<-read.table (file=paste0(args$dir,"/","plink.hwe"), header=TRUE)
pdf(paste0(args$dir,"/","histhwe.pdf"))
hist(hwe[,9],main="Histogram HWE")
dev.off()

hwe_zoom<-read.table (file=paste0(args$dir,"/","plinkzoomhwe.hwe"), header=TRUE)
pdf(paste0(args$dir,"/","histhwe_below_theshold.pdf"))
hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only")
dev.off()
