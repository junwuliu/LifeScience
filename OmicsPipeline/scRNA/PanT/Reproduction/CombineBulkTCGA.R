
source("/public/data/backup/scRNA_Anno_CancerCell/ConfirmResult/DE_Result/Cancer_vs_NormalEpi/TCGA_TPM.DE.R")
for (i in 1:dim(Pandiff)[1]){
	if (Pandiff[i,"common"] %in% c("Common","PanSpecific","MutiSpecific")){
		Pandiff[i,'TCGA_TumorHigh'] = TCGATumorHigh(gene=rownames(Pandiff[i,]),outdir="./")
	}
}
