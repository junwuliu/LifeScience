.libPaths('')
suppressMessages(library('Seurat'))
suppressMessages(library("optparse"))
suppressMessages(library(stringr))
option_list = list(make_option("--datadir", type = "character", default = NULL, help = "scRNA Rdata"),
                   make_option("--outdir", type = "character", default = NULL, help = ""))

args <- parse_args(OptionParser(option_list=option_list))
LoadToEnvironment <- function(RData){data = get(load(RData)[1]);return(data)}

#id = strsplit(basename(args$gdata),split='.',fixed=TRUE)[[1]][1]
outdir=args$outdir

metainfo = read.table(gzfile('/public/data/scDT/SingleCellDownload/PanT/RawData/GSE156728_metadata.txt.gz'),header=T,sep='\t',row.names=1,check.names=F)
count.files = list.files(path=args$datadir,pattern='counts.txt.gz',full.names=TRUE)
print (count.files)

for (i in 1:length(count.files)){
	ids = strsplit(basename(count.files[i]),split='.',fixed=TRUE)[[1]]
	id = paste(ids[1],ids[2],sep=".")
	print (count.files[i])
	print (id)
	count.data = read.table(count.files[i],header=T,sep='\t',row.names=1)
	metainfo_sub = metainfo[colnames(count.data),]
	mca = CreateSeuratObject(count.data, project = id, min.cells = 0, min.features = 0,meta.data = metainfo_sub)
	mca[["percent.mt"]] <- PercentageFeatureSet(mca, pattern="^MT-")
	mca$donor = mca$patient
	save(mca,file=paste(args$datadir,"/",id,".T.Rdata",sep=""))
}
