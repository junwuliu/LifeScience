library("optparse")
option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains"))
args <- parse_args(OptionParser(option_list=option_list))

het <- read.table(paste0(args$dir,"/","R_check.het"), head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, paste0(args$dir,"/","fail-het-qc.txt"), row.names=FALSE)
