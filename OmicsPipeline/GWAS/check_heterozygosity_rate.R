library("optparse")
option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains"))
args <- parse_args(OptionParser(option_list=option_list))

het <- read.table(paste0(args$dir,"/","R_check.het"), head=TRUE)
pdf(paste0(args$dir,"/","heterozygosity.pdf"))
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()
