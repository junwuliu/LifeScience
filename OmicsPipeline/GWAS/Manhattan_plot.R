#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
.libPaths('')
library("optparse")
library('CMplot')

option_list = list(make_option("--input", type = "character", default = NULL, help  = "gwas result"),
                   make_option("--select_name", type = "character", default = "P", help  = "select col name for CMplot"),
                   make_option("--project", type = "character", default = "test", help  = "project id"),
                   make_option("--threshold", type = "numeric", default = "2.89e-08", help  = "threshold"),
                   make_option("--output", type = "character", default = NULL, help  = "output png"))
# select_name order should be CHR/SNP/BP/P
args <- parse_args(OptionParser(option_list=option_list))

library("qqman")
library('tidyverse')

list_col = unlist(strsplit(args$select_name,','))

results_log <- read.table(args$input, head=TRUE,na.strings=NA,sep='\t')
print (dim(results_log))
results_log = results_log %>% drop_na() ## 去除带NA的行
print (dim(results_log))

#png(args$output)
#manhattan(results_log,chr="CHR",bp="BP",p=args$thesh,snp="SNP", main = "Manhattan plot: logistic")
#dev.off()

#png(args$output)

head(results_log)
results_log2 = results_log[,list_col]

#
#if plot.type="d", SNP density will be plotted; if
#          plot.type="c", only circle-Manhattan plot will be plotted; if
#          plot.type="m",only Manhattan plot will be plotted; if
#          plot.type="q",only Q-Q plot will be plotted; if
#          plot.type=c("m","q"), Both Manhattan and Q-Q plots will be plotted.


#CMplot(results_log2,plot.type=c("m","q"),LOG10=TRUE,threshold=args$threshold,file="jpg",dpi=72,file.output=TRUE,verbose=TRUE,width=18,height=8,signal.col=NULL,file.name=args$project)
CMplot(results_log2,plot.type=c("m","q"),LOG10=TRUE,threshold=args$threshold,file="jpg",dpi=72,file.output=TRUE,verbose=TRUE,width=8,height=8,signal.col=NULL,file.name=args$project)
#CMplot(results_log,plot.type="m",LOG10=TRUE,threshold=2.89e-08,file="png",dpi=500,file.output=TRUE,verbose=TRUE,width=18,height=8,signal.col=NULL,file.name=args$output)
#dev.off()


#results_as <- read.table("assoc_results.assoc", head=TRUE)
#jpeg("assoc_manhattan.jpeg")
#manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc")
#dev.off()  


