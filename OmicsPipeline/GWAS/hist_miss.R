library("optparse")

option_list = list(make_option("--dir", type = "character", default = NULL, help  = "work dir contains plink.imiss plink lmiss"))


args <- parse_args(OptionParser(option_list=option_list))

indmiss<-read.table(file=paste0(args$dir,"/","plink.imiss"), header=TRUE)
snpmiss<-read.table(file=paste0(args$dir,"/","plink.lmiss"), header=TRUE)
# read data into R 

pdf(paste0(args$dir,"/","histimiss.pdf")) #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness") #selects column 6, names header of file

pdf(paste0(args$dir,"/","histlmiss.pdf"))
hist(snpmiss[,5],main="Histogram SNP missingness")  
dev.off() # shuts down the current device
