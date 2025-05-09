.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
library('stats') ## aov
library('dplyr')
library("harmony")
option_list = list(make_option("--indir", type = "character", default = NULL, help = "metacluster csv dir"),
                  make_option("--type", type = "character", default = NULL, help = "scRNA Rdata"),
                  make_option("--Fscore", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

indir = args$indir
type= args$type

informative_datasets = list.files("*informatvieGene.csv", path=indir,full.names = TRUE)
commonGene_datasets = list.files("*CommonGene.csv", path=indir,full.names = TRUE)

total_informative_matrix = data.frame()
obj_list = list()

for (i in 1:length(commonGene_datasets)){
	id = strsplit(basename(commonGene_datasets[i]),split='.',fixed=TRUE)[[1]][1]
	print (id)
	sub_commonGene_matrix = read.table(commonGene_datasets[i],header=T,sep=",",row.names=1)
	sub_informative_matrix = read.table(informative_datasets[i],header=T,sep=",")
	mca = CreateSeuratObject(sub_commonGene_matrix, project = id, min.cells = 0, min.features = 0)
	#mca$minicluster = paste(id,"metacluster",mca$seurat_clusters,sep="_")
	mca$datasource = id
	if (id %in% c("CRC_GSE108989","GSM4743237_CHOL_SS2","HCC_GSE140228_Smartseq2","HCC_GSE98638","HNSCC_GSE103322","LC_GSE99254","MELA_GSE115978","MELA_GSE120575")){
		mca$library = "SmartSeq2"
	}else{
		mca$library = "droplet"
	}
	mca@assays$RNA@data =  as(as.matrix(sub_commonGene_matrix),"CsparseMatrix")
	print (dim(mca))
	obj_list[[i]] = mca
	if (i == 1){
		total_informative_matrix = sub_informative_matrix
	}else{
		total_informative_matrix = merge(total_informative_matrix,sub_informative_matrix,by='gene')
	}
}
list_y = obj_list[-1]
mca<-merge(obj_list[[1]],y=list_y)

rownames(total_informative_matrix) = total_informative_matrix$gene
total_informative_matrix = dplyr::select(total_informative_matrix,-gene)
mca@assays$RNA@var.features = rownames(total_informative_matrix)
mca@assays$RNA@scale.data = as.array(as.matrix(total_informative_matrix))
#object <- FindVariableFeatures(object, selection.method = select.method, nfeatures = nfeatures)
mca  <- RunPCA(mca,  npcs = 20)
#mca <- RunHarmony(mca,group.by.vars="library",dims.use = 1:15,max.iter.harmony = 10)
mca <- RunHarmony(mca,group.by.vars="datasource",dims.use = 1:15,max.iter.harmony = 10)
mca <- RunUMAP(mca,reduction = "harmony", dims=1:15)
mca <- FindNeighbors(mca,reduction = "harmony",  dims=1:15)
save(mca,file=paste("Pan",type,".T.Rdata",sep=""))

PanMarker_CD4 = list(Tn=c("LEF1","CCR7","SELL","TCF7","ADSL","IL7R"),Tcm=c("GPR183","ANXA1"),Trm=c("ZNF683","ITGAE","CXCR6","CD69"),Tem=c("GZMK","CXCR5","CCL3L3","IFNG"),Treg=c("FOXP3","RTKN2","IL2RA","CTLA4","TIGIT","TNFRSF9","S1PR1"),Temra=c("PRF1","GNLY","GZMB","GZMH","FGFBP2","KLRG1"),ISG=c("ISG15","IFIT3"),Other=c("NME1","CCR4","CCL5","CXCL13","AREG","CXCR5","CREM","IL21","IL26"))
PanMarker_CD8 = list(Tex=c("LAG3","PDCD1","HAVCR2",'TOX','TOX2','ENTPD1','LAYN','PRDM1',"ID2"),Tem=c("GZMK","CXCR5","CCL3L3","IFNG"),Tn=c("LEF1","CCR7","SELL","TCF7"),Tcm=c("GPR183","ANXA1"),Trm=c("ZNF683","ITGAE","CXCR6","CD69"),MAIT=c("SLC4A10", "RAV1-2","KLRB1"),Temra=c("PRF1","GNLY","GZMB","GZMH","FGFBP2","KLRG1"),NK_like=c("KIR3DL1"),ISG=c("ISG15","IFIT3"),Tc17=c("RORC","KLRB1","IL23R","CCR6"),Other=c("NME1","CXCL13","FOXP3","EOMES","OXPHOS","TXK"))


