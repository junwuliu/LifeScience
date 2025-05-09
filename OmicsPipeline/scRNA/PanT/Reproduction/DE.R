source('/public/scRNA_pipline/Module/DEGene.R')
 diffGene = DEGene(mca,metaname='subC',case="Terminal_Tex",control="TCF7_Tex")
diffGene = diffGene[order(-diffGene$log2FC),]
sigGene = subset(diffGene,abs(log2FC)>=0.25 & Adj_pval<0.05)


source("/public/scRNA_pipline/Module/UniprotSubcellularLocation.R")
for (i in 1:dim(sigGene)[1]){
	if (is.na(sigGene[i,"Location"])){
		sigGene[i,"Location"] = locations(sigGene[i,"gene"])
	}
}
