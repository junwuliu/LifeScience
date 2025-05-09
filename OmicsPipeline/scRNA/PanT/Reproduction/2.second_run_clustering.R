.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
library('stats') ## aov
library('dplyr')
option_list = list(make_option("--data", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--library", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

#outdir=args$outdir

source('/public/data/scDT/SingleCellDownload/PanT/Script/Nomalization.R')
mca = LoadToEnvironment(args$data)
base_name = basename(args$data)
id = strsplit(basename(args$data),split='.',fixed=TRUE)[[1]][1]

immunoglobulin_gene = read.table("/public/data/scDT/SingleCellDownload/PanT/Script/immunoglobulin.genelist",header=F,sep="\t")$V1
TCR_gene = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/TCR.genelist',header=F,sep='\t')$V1
proliferation_gene = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/proliferation.genelist',header=F,sep='\t')$V1
DIG_genes  = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/Dissocation_gene.list',header=F,skip=1,sep='\t')$V1
DIG = list()
DIG[[1]] = DIG_genes
DIG_genes_new = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/Dissocation_gene.FromDatasets.list.merge',header=F,sep='\t')$V1
DIG[[2]] = DIG_genes_new

mca = NormalSelf(mca,library=args$library)
if (args$library == 'Smart'){
	data_matrix = log2(mca@assays$RNA@data+1)
	DefaultAssay(mca) = 'RNA'
	mca@assays$RNA@data = as(data_matrix,"CsparseMatrix")
	mca = FindVariableFeatures(mca, selection.method = 'vst', nfeatures = 6000)
}else{
	data_matrix = log2(mca@assays$RNA@data+1)
	DefaultAssay(mca) = 'RNA'
	mca@assays$RNA@data = as(data_matrix,"CsparseMatrix")
	mca = FindVariableFeatures(mca, selection.method = 'vst', nfeatures = 6000)
}
variableGene = mca@assays$RNA@meta.features
print (dim(variableGene))
#variableGene = sort(variableGene[,'vst.variance.standardized'],decreasing=T)
variableGene = variableGene[order(variableGene$vst.variance.standardized,decreasing=T),]

head(variableGene)
RP_gene = grep(pattern = "^RP([0-9]+-|[LS])",rownames(variableGene),perl = TRUE,value=TRUE)
variableGene = variableGene[!(rownames(variableGene) %in% c(immunoglobulin_gene,TCR_gene,proliferation_gene,RP_gene,"MALAT1",DIG_genes,DIG_genes_new)),]
print (dim(variableGene))
head(variableGene)
variableGene = variableGene %>% top_n("vst.variance.standardized", n = 1500)

cc.genes = LoadToEnvironment('/public/data/scDT/SingleCellDownload/PanT/Script/cc.genes.updated.2019.rda')
mca <- CellCycleScoring(mca, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
mca <- AddModuleScore(mca,features=DIG,name = "DIG.Score")
head(mca[])
#for (i in 0:(length(unique(mca$seurat_clusters))-1)){tmp = mca[,mca$seurat_clusters == i];print (median(tmp$DIG.Score1))}
if (length(unique(mca$donor)) > 1){
	mca <- ScaleData(mca,features = rownames(variableGene),vars.to.regress = c('donor',"S.Score", "G2M.Score","nCount_RNA",'percent.mt','DIG.Score1',"DIG.Score2"))
}else{
	mca <- ScaleData(mca,features = rownames(variableGene),vars.to.regress = c("S.Score", "G2M.Score","nCount_RNA",'percent.mt','DIG.Score1',"DIG.Score2"))
}

mca <- RunPCA(mca, npcs = 15)
mca = RunUMAP(mca, dims=1:15)
mca <- FindNeighbors(mca,dims=1:15)
if (dim(mca)[2] > 50){
	mca = FindClusters(mca,resolution = 2)
}else{
	mca = FindClusters(mca,resolution = 1)
}

#load(paste(args$outdir,"/",base_name,sep=''))
source("/public/scRNA_pipline/Module/marker.R")
if (args$library == 'Smart'){
	logFCfilter = 1
}else{
	logFCfilter = 0.25
}
#top_marker_list = AllTopMarker(mca,idents="seurat_clusters",resolution=2,outdir=args$outdir,topnum=30,prefix=base_name,logFCfilter=logFCfilter,adjPvalFilter=0.01,min.pct=0.25,rmMT=FALSE,assays="RNA")


## ANOVA test
print ("start to run ANOVA test")
data_matrix = mca@assays$RNA@data
metainfo = mca@meta.data
total_gene_num = dim(data_matrix)[1]

gene_Fscore = data.frame()
for (gene in rownames(data_matrix)){
	gene_matrix = data.frame(data_matrix[gene,])
	names(gene_matrix) = "data"
	gene_matrix$seurat_clusters = metainfo[match(rownames(gene_matrix),rownames(metainfo)),'seurat_clusters']
	aov_result = summary(aov(data~seurat_clusters,data = gene_matrix))
	F_score = aov_result[[1]]["F value"][[1]][1]
	gene_F = data.frame(gene,F_score)
	names(gene_F) = c('gene',"F_score")
	gene_Fscore = rbind(gene_Fscore,gene_F)
}
gene_Fscore = gene_Fscore[order(gene_Fscore$F_score,decreasing=T),]
for (i in 1:dim(gene_Fscore)[1]){
	gene_Fscore[i,"PercentileRank"] = round(i / total_gene_num,10)
}

write.table(gene_Fscore,file=paste(args$outdir,"/",id,".GeneFscore.xls",sep=''),sep="\t",quote=F,row.names=F,col.names=T)

save(mca,file=paste(args$outdir,"/",base_name,sep=''))
