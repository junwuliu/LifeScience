Limma <- function(object=NULL,metaname=NULL,case=NULL,control=NULL,maxNum=0,selectGene=NULL){
	if (is.null(object) || is.null(metaname) || is.null(case) || is.null(control)){
		writeLines ("Exam: diffGene = DEGene(mca,metaname='Gender',case='case',control='control')
        Parameter:
            object: Seurat object name
            metaname: which colname of metadata contains case and control
            case: case name
            control: control name
        ")
    stop('Please give a Seurat Object name')
	}
	library('limma')
	library("dplyr")
	library("ggplot2")
	library("ggrepel")
	library('tools')
	library('metaMA')
	pseudocount.use =1
	Group1 <- rownames(object[,object@meta.data[,metaname] %in% case]@meta.data)
	case_cellnum = length(Group1)
	print (paste("Group ",case," meta-clusters cells are ",case_cellnum,sep=""))
	Group2 <- rownames(object[,object@meta.data[,metaname] %in% control]@meta.data)
	control_cellnum = length(Group2)
	if (maxNum>0 & control_cellnum >= maxNum){
		Group2 = sample(Group2,maxNum)
		print (paste("Group ",control," meta-clusters cells are ",control_cellnum," ,subset to",maxNum,sep=""))
	}else{
		print (paste("Group ",control," meta-clusters cells are ",control_cellnum,sep=""))
	}
	Expression_1 <- as.matrix(object@assays$RNA@data[,Group1])
	Expression_2 <- as.matrix(object@assays$RNA@data[,Group2])
	if (!is.null(selectGene)){
		Expression_1 = Expression_1[selectGene,]
		Expression_2 = Expression_2[selectGene,]
	}
	#Expression_1 = 2^Expression_1 - 1
	#Expression_2 = 2^Expression_2 - 1
	print (Expression_1[1:2,1:2])
	## limma
	MergeMatrix = cbind(Expression_1,Expression_2)
	list_group = c(rep(case,dim(Expression_1)[2]),rep(control,dim(Expression_2)[2])) %>% factor(.,levels=c(case,control), ordered = F )
	list_group <- model.matrix(~factor(list_group)+0)
	colnames(list_group) = c("cases","controls")
	df.fit <- lmFit(MergeMatrix, list_group)
	df.matrix <- makeContrasts(cases - controls , levels = list_group)
	fit <- contrasts.fit(df.fit, df.matrix)
	fit <- eBayes(fit,robust = TRUE)
	#tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
  ## 自带 BH correction, 无需再做
  ## effect size table
  result = effectsize(fit$t,dim(MergeMatrix)[2],(fit$df.prior+fit$df.residual))
  #result = effectsize(fit$t,(fit$df.prior+fit$df.residual),dim(MergeMatrix)[2])
  result <- data.frame(result)
  result = result[order(result$dprime,decreasing=T),]
  return(result)
}

