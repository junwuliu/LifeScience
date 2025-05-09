library("optparse")
option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains"))
args <- parse_args(OptionParser(option_list=option_list))

gender <- read.table(paste0(args$dir,"/","plink.sexcheck"), header=T,as.is=T)

pdf(paste0(args$dir,"/","Gender_check.pdf"))
hist(gender[,6],main="Gender", xlab="F")
dev.off()

pdf(paste0(args$dir,"/","Men_check.pdf"))
male=subset(gender, gender$PEDSEX==1)
hist(male[,6],main="Men",xlab="F")
dev.off()

pdf(paste0(args$dir,"/","Women_check.pdf"))
female=subset(gender, gender$PEDSEX==2)
hist(female[,6],main="Women",xlab="F")
dev.off()

