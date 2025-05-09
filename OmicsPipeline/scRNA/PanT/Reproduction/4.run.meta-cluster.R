.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
library('stats') ## aov
library('dplyr')
option_list = list(make_option("--data", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--Fscore", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--library", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

outdir=args$outdir
mca = LoadToEnvironment(args$data)
id = strsplit(basename(args$data),split='.',fixed=TRUE)[[1]][1]

## 先获得根据所有基因scale的矩阵，然后得到每一个minicluster的平均值矩阵做为PanData的data矩阵
# 获得根据informative gene的scale的矩阵，并作为后续PanData的scale矩阵用于分群
# Thus, the original gene by cell expression matrix was converted to the gene by mini-cluster expression matrix

totalscore = read.table(args$Fscore,header=T,sep="\t",row.names=1)
commonGenes = rownames(totalscore)
totalscore = totalscore[order(totalscore$median,decreasing=F),]

informative_genes = rownames(head(totalscore,1500))

### get meta clusters
print ("start to run meta clusters")
if ("Donor" %in% colnames(mca@meta.data)){
	mca$donor = mca$Donor
}

if (length(unique(mca$donor)) > 1){
	mca <- ScaleData(mca,features = commonGenes, vars.to.regress = c('donor',"S.Score", "G2M.Score","nCount_RNA",'percent.mt','DIG.Score1',"DIG.Score2"))
	scale_matrix = mca@assays$RNA@scale.data
	mca <- ScaleData(mca,features = informative_genes, vars.to.regress = c('donor',"S.Score", "G2M.Score","nCount_RNA",'percent.mt','DIG.Score1',"DIG.Score2"))
}else{
	mca <- ScaleData(mca,features = commonGenes, vars.to.regress = c("S.Score", "G2M.Score","nCount_RNA",'percent.mt','DIG.Score1',"DIG.Score2"))
	scale_matrix = mca@assays$RNA@scale.data
	mca <- ScaleData(mca,features = informative_genes, vars.to.regress = c("S.Score", "G2M.Score","nCount_RNA",'percent.mt','DIG.Score1',"DIG.Score2"))
}
# 
#mca@assays$RNA@data = as.array(as.matrix(scale_matrix))

mca <- RunPCA(mca, npcs = 15)
mca = RunUMAP(mca, dims=1:15)
mca <- FindNeighbors(mca,dims=1:15,k.param=10)
if (dim(mca)[2] < 50){
	mca = FindClusters(mca,resolution = 1)
}else if (dim(mca)[2] < 500){
	mca = FindClusters(mca,resolution = 25)
}else{
	mca = FindClusters(mca,resolution = 50)
}
save(mca,file=paste(args$outdir,"/",id,".miniclusters.Rdata",sep=''))
#break
meta_cluster = data.frame(table(mca$seurat_clusters))
names(meta_cluster) = c("clusters","cellnum")

write.table(meta_cluster,file=paste(args$outdir,"/",id,".meta_clusters.xls",sep=''),sep="\t",quote=F,row.names=F,col.names=T)
### 
print ("start to creat meta-cluster matrix")
#

informative_matrix = mca@assays$RNA@scale.data
total_mini_cluster = data.frame()
total_informative_cluster = data.frame()
for (i in 1:length(unique(mca$seurat_clusters))){
	mini_cluster = mca[,mca$seurat_clusters == i-1]
	mini_cluster_cell = rownames(mini_cluster@meta.data)
	mini_cluster_matrix = data.frame(scale_matrix[,mini_cluster_cell])
	mini_cluster_informative = data.frame(informative_matrix[,mini_cluster_cell])
	cluster_name = paste(id,".minicluster_",i-1,sep="")
	mini_cluster_informative[,cluster_name] = apply(mini_cluster_informative,1,mean)
	mini_cluster_informative$gene = rownames(mini_cluster_informative)
	mini_cluster_matrix[,cluster_name] = apply(mini_cluster_matrix,1,mean)
	mini_cluster_matrix$gene = rownames(mini_cluster_matrix)
	if (i == 1){
		total_informative_cluster = mini_cluster_informative[,c("gene",cluster_name)]
		total_mini_cluster = mini_cluster_matrix[,c("gene",cluster_name)]
		
	}else{
		total_mini_cluster = merge(total_mini_cluster,mini_cluster_matrix[,c("gene",cluster_name)],by="gene")
		total_informative_cluster = merge(total_informative_cluster,mini_cluster_informative[,c("gene",cluster_name)],by="gene")
	}
}
#rownames(total_mini_cluster) = total_mini_cluster$gene
#total_mini_cluster = select(total_mini_cluster,-gene)
write.table(total_mini_cluster,file=paste(args$outdir,"/",id,".meta_clusters.CommonGene.csv",sep=''),sep=",",quote=F,row.names=F,col.names=T)
write.table(total_informative_cluster,file=paste(args$outdir,"/",id,".meta_clusters.informatvieGene.csv",sep=''),sep=",",quote=F,row.names=F,col.names=T)
print ("Proecess done")

