.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
option_list = list(make_option("--data", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

#outdir=args$outdir

mca = LoadToEnvironment(args$data)
id = strsplit(basename(args$data),split='.',fixed=TRUE)[[1]][1]

CD3_pos = subset(mca,CD3D >0 | CD3G > 0)
print (paste("total T is",dim(mca)[2],";CD3 Pos is",dim(CD3_pos)[2]),sep="")

CD4_pos = subset(CD3_pos, CD4 >0 & (CD8A == 0 & CD8B ==0))
CD8_pos = subset(CD3_pos, CD4 ==0 & (CD8A >0 | CD8B >0))
print (paste("CD4 pos is",dim(CD4_pos)[2],";CD8 Pos is",dim(CD8_pos)[2]),sep="")
print (paste(id,dim(CD4_pos)[2],dim(CD8_pos)[2]),sep=" ")

save(CD4_pos,file=paste(args$outdir,"/",id,".CD4.T.Rdata",sep=""))
save(CD8_pos,file=paste(args$outdir,"/",id,".CD8.T.Rdata",sep=""))

