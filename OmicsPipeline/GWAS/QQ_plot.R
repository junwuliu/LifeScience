#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
.libPaths('')
library("optparse")

option_list = list(make_option("--input", type = "character", default = NULL, help  = "gwas result"),
                  make_option("--output", type = "character", default = NULL, help  = "QQ-plot png"))
args <- parse_args(OptionParser(option_list=option_list))

library("qqman")
library('CMplot')

results_log <- read.table(args$input, head=TRUE,sep='\t',na.strings='')
png(args$output)
qq(results_log$FDR_BH, main = "Q-Q plot of GWAS p-values : log")
dev.off()

#results_as <- read.table("assoc_results.assoc", head=TRUE)
#jpeg("QQ-Plot_assoc.jpeg")
#qq(results_as$P, main = "Q-Q plot of GWAS p-values : log")
#dev.off()

