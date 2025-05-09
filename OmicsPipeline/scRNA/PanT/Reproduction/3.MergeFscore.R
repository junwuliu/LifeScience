.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
library('stats') ## aov
library('dplyr')
option_list = list(make_option("--indir", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

outdir=args$outdir
indir=args$indir

totalscore = data.frame()
Fscorelist = list.files('*GeneFscore.xls',path=indir,full.names =TRUE)
print (length(Fscorelist))
for (i in 1:length(Fscorelist)){
	score = read.table(Fscorelist[i],header=T,sep="\t")
	score = score[,c("gene","PercentileRank")]
	id = strsplit(basename(Fscorelist[i]),split='.',fixed=TRUE)[[1]][1]
	names(score) = c("gene",id)
	if (i == 1){
		totalscore = score
	}else{
		totalscore = merge(totalscore,score)
	}
}
immunoglobulin_gene = read.table("/public/data/scDT/SingleCellDownload/PanT/Script/immunoglobulin.genelist",header=F,sep="\t")$V1
TCR_gene = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/TCR.genelist',header=F,sep='\t')$V1
proliferation_gene = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/proliferation.genelist',header=F,sep='\t')$V1
DIG_genes  = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/Dissocation_gene.list',header=F,skip=1,sep='\t')$V1
DIG_genes_new = read.table('/public/data/scDT/SingleCellDownload/PanT/Script/Dissocation_gene.FromDatasets.list',header=F,sep='\t')$V1
RP_gene = grep(pattern = "^RP([0-9]+-|[LS])",totalscore$gene,perl = TRUE,value=TRUE)

totalscore = totalscore[!totalscore$gene %in% c(immunoglobulin_gene,TCR_gene,proliferation_gene,DIG_genes,DIG_genes_new,RP_gene),]
rownames(totalscore) = totalscore$gene
totalscore = select(totalscore,-gene)
totalscore$median = apply(totalscore,1,median)
totalscore = totalscore[order(totalscore$median,decreasing=F),]
write.table(totalscore,file=paste(outdir,"MergeAllDatasets.Fscore.xls",sep="/"),sep="\t",row.names=T,col.names=T,quote=F)
#totalscore = totalscore[order(totalscore$median,decreasing=F),] %>% top_n("median",1500)
