NormalSelf<-function(object=NULL,library=NULL){
	if(is.null(object) || is.null(library)){
        writeLines ("回溯亚群的meta列内容返回到到原来的object,根据cell_id一致性
        Exam: raw=TraceBack(raw=mca,sub=subT,meta='subC'))
        Parameter:
             sub: sub seurat object name
             meta: colnames of sub object and raw object
        ")
        stop('please give raw object and sub object')
	}
	library('Seurat')
	if (library == 'Smart'){
		library(BiocParallel)
		library(scran)
		sce <- as.SingleCellExperiment(object)
		clusters <- quickCluster(sce)
		sce <- computeSumFactors(sce, clusters=clusters)
		summary(sizeFactors(sce))
		#counts(sce) <- as(counts(sce, withDimnames = FALSE), "dgCMatrix")
		object <- as.Seurat(sce, counts = "counts", data = "logcounts")
		return (object)
	}else{
		#object <- SCTransform(object, verbose = FALSE)
		print ("using NormalizeData")
		object <- NormalizeData(object, normalization.method = 'LogNormalize', scale.factor = 10000)
	}
	return (object)
}
