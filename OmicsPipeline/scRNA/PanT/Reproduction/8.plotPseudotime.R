.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
library('dplyr')
library('reshape2')
library('ggplot2')
option_list = list(make_option("--data", type = "character", default = NULL, help = "monocle rdata"),
                  make_option("--Rdata", type = "character", default = NULL, help = "panT rdata"),
                  make_option("--gene", type = "character", default = NULL, help = "gene name"),
                  make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

cds = LoadToEnvironment(args$data)
mca =  LoadToEnvironment(args$Rdata)
gene = args$gene

pseudo = data.frame(cds@principal_graph_aux@listData$UMAP$pseudotime)
names(pseudo) = c("pseudotime")
metadata = mca@meta.data
pseudo$subC = metadata[match(rownames(pseudo),rownames(metadata)),"newSubC"]

data_matrix = mca@assays$RNA@data
gene_matrix = data.frame(data_matrix[gene,])
names(gene_matrix) = c("RelativeExp")
pseudo$exp = gene_matrix[match(rownames(pseudo),rownames(gene_matrix)),"RelativeExp"]


# select = object[,object$seurat_clusters %in% c(0,1,2,3,4,8)]
#new_pseudo = pseudo[rownames(object@meta.data),]
#new_pseudo$exp = gene_matrix[match(rownames(new_pseudo),rownames(gene_matrix)),"RelativeExp"]
png(file=paste(args$outdir,"/",gene,".pseudotime.expr.png",sep=""),width = 600, height = 400)
## 分组点图
#ggplot(data=pseudo, aes(x=pseudotime, y=exp))+ geom_point(aes(color=subC),alpha=1)+ geom_smooth(method = "loess") + ylim(-0.5, 1) + ggtitle(gene)
## loess 拟合图
ggplot(data=pseudo, aes(x=pseudotime, y=exp))+ geom_smooth(method = "loess") + ylim(-0.5, 1) + ggtitle(gene)
dev.off()

