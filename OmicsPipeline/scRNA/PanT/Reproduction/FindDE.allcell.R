.libPaths('')
load("/public/data/scDT/SingleCellDownload/PanT/1_selectData/outdir/step2_cluster/Exclude_LC_GSE117570/PanCD8.T.Rdata")
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

totalscore = read.table("/public/data/scDT/SingleCellDownload/PanT/1_selectData/outdir/step2_cluster/Exclude_LC_GSE117570/MergeAllDatasets.Fscore.xls",header=T,sep="\t",row.names=1)
commonGenes = rownames(totalscore)

TerminalTex = mca[,mca$newSubC %in% c("Terminal_Tex")]
TCF7Tex = mca[,mca$newSubC %in% c("TCF7_Tex")]

TerminalTex_metacluster = rownames(TerminalTex@meta.data)
TCF7Tex_metacluster = rownames(TCF7Tex@meta.data)


Terminal_Tex_datamatrix = matrix()
TCF7_Tex_datamatrix = matrix()


allTerminalTex_datasource = unique(TerminalTex$datasource)
allTCF7Tex_datasource = unique(TCF7Tex$datasource)
for (i in 1:length(allTerminalTex_datasource)){
	id = allTerminalTex_datasource[i]
	print (id)
	minicluters = rownames(TerminalTex[,TerminalTex$datasource == id]@meta.data)
	subData = LoadToEnvironment(paste("/public/data/scDT/SingleCellDownload/PanT/1_selectData/outdir/step2_cluster/Exclude_LC_GSE117570/",id,".miniclusters.Rdata",sep=""))
	id = gsub("-",".",id)
	subData$miniclusters = paste(id,".minicluster","_",subData$seurat_clusters,sep="")
	subcluster = subData[,subData$miniclusters %in% minicluters]
	datamatrix = subcluster@assays$RNA@data
	datamatrix = datamatrix[commonGenes,]
	#subcluster@assays$RNA@data = datamatrix
	subcluster$PanCellType = "TerminalTex"
	if (i == 1){
		Terminal_Tex_datamatrix = datamatrix
	}else{
		Terminal_Tex_datamatrix = cbind(Terminal_Tex_datamatrix,datamatrix)
	}
}

print ("start to find TCF7")
for (i in 1:length(allTCF7Tex_datasource)){
	id = allTCF7Tex_datasource[i]
	print (id)
	minicluters = rownames(TCF7Tex[,TCF7Tex$datasource == id]@meta.data)
	subData = LoadToEnvironment(paste("/public/data/scDT/SingleCellDownload/PanT/1_selectData/outdir/step2_cluster/Exclude_LC_GSE117570/",id,".miniclusters.Rdata",sep=""))
	id = gsub("-",".",id)
	subData$miniclusters = paste(id,".minicluster","_",subData$seurat_clusters,sep="")
	subcluster = subData[,subData$miniclusters %in% minicluters]
	subcluster$PanCellType = "TCF7Tex"
	datamatrix = subcluster@assays$RNA@data
    datamatrix = datamatrix[commonGenes,]
	if (i == 1){
        TCF7_Tex_datamatrix = datamatrix
    }else{
        TCF7_Tex_datamatrix = cbind(TCF7_Tex_datamatrix,datamatrix)
    }
}

Expression_1 = 2^Terminal_Tex_datamatrix - 1
Expression_2 = 2^TCF7_Tex_datamatrix - 1

pseudocount.use = 1
mean_c1 <- as.data.frame(rowMeans(expm1(as.matrix(Expression_1))))
colnames(mean_c1) <- "mean_c1"
mean_c2 <- as.data.frame(rowMeans(expm1(as.matrix(Expression_2))))
colnames(mean_c2) <- "mean_c2"

log2fc <- data.frame(log2fc = log2(mean_c1$mean_c1 + pseudocount.use) - log2(mean_c2$mean_c2 + pseudocount.use))
rownames(log2fc) <- rownames(mean_c1)
log2fc$gene <- rownames(log2fc)

Group1 = colnames(Expression_1)
Group2 = colnames(Expression_2)
group.info <- data.frame(row.names = c(Group1, Group2))
group.info[, "group"] <- factor(x = group.info[, "group"])
data.use = cbind(Terminal_Tex_datamatrix,TCF7)
data.use = 2^data.use - 1
p_val <- sapply(
	X = 1:nrow(x = data.use),FUN = function(x) {
    return(wilcox.test(data.use[x, ] ~ group.info[, "group"])$p.value)
})

