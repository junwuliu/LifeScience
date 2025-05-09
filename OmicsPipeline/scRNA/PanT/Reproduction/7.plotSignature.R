.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
library('dplyr')
library('reshape2')
library('ggplot2')
option_list = list(make_option("--indir", type = "character", default = NULL, help = "dir include all sub celltype signature xls"),
                  make_option("--gene", type = "character", default = NULL, help = "gene name"),
                  make_option("--type", type = "character", default = NULL, help = "CD4 or CD8"),
                  make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

genes = unlist(strsplit(args$gene,split=','))
indir=args$indir
outdir=args$outdir
if (args$type == "CD8"){
	orders = c("Tn","IL7R+ Tm","ZNF683+CXCR6- Trm","ZNF683+CXCR6+ Trm","Temra","KIR+EOMES+ NK-like","KIR+TXK+ NK-like","Early Tem","Tem","GZMK+ Tex","OXPHOS_Tex","Terminal_Tex","TCF7_Tex","FOXP3_Tex","ISG+ T","Tc17","NME1_T","Un")
}else{
	orders = c("Tn","ADSL+ Tn","IL7R- Tn","TNF+ T","ISG+ T","Tm","AREG+ Tm","TIMP1 Tm","CREM+ Tm","CAPG+ Tm","GZMK+ Tem","Temra","Th17","Tfh","Tfh_Th","TNFRSF9- Treg","TNFRSF9+ Treg","Un")
}
max_es = 0.51
min_es = -0.26

sigs = list.files("signatureGene.xls$",path=indir,full.names=T)

all_data = data.frame()
for (i in 1:length(sigs)){
	es = read.table(sigs[i],header=T,sep="\t",row.names=1)
	id = strsplit(basename(sigs[i]),split='.',fixed=TRUE)[[1]][1]
	es = es[,c("gene","avg_ES")]
	names(es) = c('gene',id)
	if (i==1){
		all_data = es
	}else{
		all_data = merge(all_data,es,by='gene')
	}
}
rownames(all_data) = all_data$gene


## ggplot
plot_data = all_data[genes,]
plot_data = reshape2::melt(plot_data)
names(plot_data) = c("gene","cellType","ES")
plot_data$gene = factor(plot_data$gene,levels=rev(genes))
plot_data$cellType = factor(plot_data$cellType,levels=orders)

plot_data$ES[plot_data$ES > max_es] = max_es
plot_data$ES[plot_data$ES < min_es] = min_es

p = ggplot(data=plot_data, aes(x=cellType, y=gene))+ geom_point(aes(color=ES,size=ES),alpha=1)+ theme_classic()+ theme(panel.grid.major = element_line(color = "gray"),panel.border = element_rect(color="black",fill=NA))+ theme(axis.text.x = element_text(size=12,angle=60,hjust=1),axis.text.y =element_text(size=12))+ xlab(NULL)+ylab(NULL)+ scale_colour_distiller(palette = "RdYlBu")

### plot TCF7 high ES
#all_data = read.table("/public/data/scDT/SingleCellDownload/PanT/1_selectData/outdir/step2_cluster/Exclude_LC_GSE117570/ClusterCompare/CD8T.MergedSignature.xls",header=T,sep="\t",row.names=1)
#TCF7Tex_high = all_data[order(all_data$TCF7_Tex,decreasing = T),]
#TCF7_mem = TCF7Tex_high[grep("membrane",TCF7Tex_high$Location,ignore.case = T),]
#genes = c(head(TCF7_mem$gene,30),"PDCD1")

#plot_data = all_data[genes,]
#plot_data = reshape2::melt(plot_data)
#names(plot_data) = c("gene",'location',"cellType","ES")
#plot_data$gene = factor(plot_data$gene,levels=rev(genes))
#plot_data$cellType = factor(plot_data$cellType,levels=orders)


