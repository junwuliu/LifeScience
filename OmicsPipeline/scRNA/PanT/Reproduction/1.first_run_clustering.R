.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
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

mca = NormalSelf(mca,library=args$library)
#countMatrix = mca@assays$RNA@counts
#if (exists(countMatrix)){
#	mca = NormalSelf(mca,library=args$library)
#}else{
#	print ('Dont have count matrix, use CPM data matrix')
#}
# 对count重新进行归一化
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
## 获得高可变基因，去除一些干扰基因
variableGene = variableGene[order(variableGene$vst.variance.standardized,decreasing=T),]

head(variableGene)
RP_gene = grep(pattern = "^RP([0-9]+-|[LS])",rownames(variableGene),perl = TRUE,value=TRUE)
variableGene = variableGene[!(rownames(variableGene) %in% c(immunoglobulin_gene,TCR_gene,proliferation_gene,RP_gene,"MALAT1")),]
print (dim(variableGene))
head(variableGene)
variableGene = variableGene %>% dplyr::top_n(vst.variance.standardized, n = 1500)

cc.genes = LoadToEnvironment('/public/data/scDT/SingleCellDownload/PanT/Script/cc.genes.updated.2019.rda')
mca <- CellCycleScoring(mca, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

if ("Donor" %in% colnames(mca@meta.data)){
	mca$donor = mca$Donor
}

if (length(unique(mca$donor)) > 1){
	mca <- ScaleData(mca,features = rownames(variableGene),vars.to.regress = c('donor',"S.Score", "G2M.Score"))
}else{
	mca <- ScaleData(mca,features = rownames(variableGene),vars.to.regress = c("S.Score", "G2M.Score"))
}

mca <- RunPCA(mca, npcs = 15)
mca = RunUMAP(mca, dims=1:15)
mca <- FindNeighbors(mca,dims=1:15)
if (dim(mca)[2] > 50){
	mca = FindClusters(mca,resolution = 2)
}else{
	mca = FindClusters(mca,resolution = 1)
}

save(mca,file=paste(args$outdir,"/",base_name,sep=''))
#load(paste(args$outdir,"/",base_name,sep=''))
source("/public/scRNA_pipline/Module/marker.R")
if (args$library == 'Smart'){
	logFCfilter = 1
}else{
	logFCfilter = 0.25
}
if (dim(mca)[2] > 50){
top_marker_list = AllTopMarker(mca,idents="seurat_clusters",resolution=2,outdir=args$outdir,topnum=30,prefix=base_name,logFCfilter=logFCfilter,adjPvalFilter=0.01,min.pct=0.25,rmMT=FALSE,assays="RNA")
}else{
top_marker_list = AllTopMarker(mca,idents="seurat_clusters",resolution=1,outdir=args$outdir,topnum=30,prefix=base_name,logFCfilter=logFCfilter,adjPvalFilter=0.01,min.pct=0.25,rmMT=FALSE,assays="RNA")
}


