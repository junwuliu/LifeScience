library("optparse")
option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains"))
args <- parse_args(OptionParser(option_list=option_list))

maf_freq <- read.table(paste0(args$dir,"/","MAF_check.frq"), header =TRUE, as.is=T)
pdf(paste0(args$dir,"/","MAF_distribution.pdf"))
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
dev.off()


