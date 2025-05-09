.libPaths('')
library("optparse")
library(ggpubr)
library(ggrepel)
library(patchwork)
library(ggplot2)
library(gggenes)
library(gridExtra)
library('extrafont')
library(tidyr)
library(dplyr)

option_list = list(make_option("--gene", type = "character", default = NULL, help = "genelist info"),
                    make_option("--gpd", type = "character", default = NULL, help = "gene pred annotation file"),
                    make_option("--gscore", type = "character", default = NULL, help = "gRNA score file"),
                    make_option("--outfile", type = "character", default = NULL, help = "out of png file"))

args <- parse_args(OptionParser(option_list=option_list))

outfile=args$outfile
gene = args$gene
gpd_file = args$gpd
gscore_file = args$gscore

gpd = read.table(gpd_file,header=T,sep='\t')
gpd$transcript_len = gpd$txEnd - gpd$txStart
gpd = gpd[order(gpd$transcript_len,decreasing=T),]
if (dim(gpd)[1] > 10){
	gpd = gpd[1:10,]
}
plot_data = list()
cds_data = list()
for (i in 1:nrow(gpd)){
	chr = gpd[i,'chr']
	genename = gpd[i,'geneName']
	transcript_id = gpd[i,'transcript_id']
	exonstarts = unlist(strsplit(gpd[i,'exonStarts'],","))
	exonends = unlist(strsplit(gpd[i,'exonEnds'],","))
	cds_s = gpd[i,'cdsStart']
	cds_e = gpd[i,'cdsEnd']
	strand = gpd[i,'strand']
	cds_data = append(cds_data,list(c(transcript_id,cds_s,cds_e,strand)))
	total_exonnum = length(exonstarts)
	direction = ifelse(strand == "+", 1,0 )
	for (j in 1:total_exonnum){
		exon_num = ifelse(strand == "+", j, total_exonnum-j+1)
		plot_data = append(plot_data,list(c(chr,'exon',as.numeric(exonstarts[j])+1,exonends[j],strand,genename,transcript_id,exon_num,direction)))
	}
}
merged_plot = data.frame(do.call(rbind, plot_data))
names(merged_plot) = c('CHROM','Type','start','end','strand','gene_id','transcript','exon_number','direction')
merged_cds = data.frame(do.call(rbind,cds_data))
names(merged_cds) = c('transcript','cds_start','cds_end','strand')
merged_plot$start = as.numeric(merged_plot$start)
merged_plot$end = as.numeric(merged_plot$end)
merged_plot$direction = as.numeric(merged_plot$direction)


min_x = min(merged_plot$start) - 500
max_x = max(merged_plot$end) + 500

### plot Structure
title = unique(gpd$geneName)
length_transcript = length(unique(gpd$transcript_id))
if (length_transcript <=12 ){
p1 = ggplot(merged_plot, aes(xmin = start, xmax = end, y = strand, fill = transcript ,label=exon_number,forward = direction)) +
    #geom_gene_arrow(arrow_body_height = grid::unit(4, "mm"),arrowhead_width = grid::unit(5, "mm"),arrowhead_height = grid::unit(5, "mm")) + geom_gene_label(aes(label=exon_number),size=8) +
	geom_gene_arrow() + geom_gene_label(aes(label=exon_number),size=8) +
    facet_grid(transcript ~ .) + labs(x=NULL,y=NULL) + ## 一个坐标轴
    scale_fill_brewer(palette = "Set3") +
    xlim(c(min_x,max_x)) +
    theme_genes()  + theme(axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank()) +
    labs(title=title) + theme(plot.title = element_text(size=16, face="bold",vjust=0)) +
    theme(legend.position = "top",legend.text=element_text(size=10),legend.key.size=unit(0.2,"inches"))+guides(color = guide_legend(ncol = NULL, byrow = TRUE) + theme(legend.spacing.x=unit(3,'cm'),legend.border.width=20))
}else{
#gradient_colors <- gradient(n = 100, colspace = "husl", space = "rgb")
color_palette <- colorRampPalette(c("blue", "red"))
gradient_colors <- color_palette(100)
p1 = ggplot(merged_plot, aes(xmin = start, xmax = end, y = strand, fill = transcript ,label=exon_number,forward = direction)) +
    #geom_gene_arrow(arrow_body_height = grid::unit(4, "mm"),arrowhead_width = grid::unit(5, "mm"),arrowhead_height = grid::unit(5, "mm")) + geom_gene_label(aes(label=exon_number),size=8) +
	geom_gene_arrow() + geom_gene_label(aes(label=exon_number),size=8) +
    facet_grid(transcript ~ .) + labs(x=NULL,y=NULL) + ## 一个坐标轴
    scale_color_gradientn(colors = gradient_colors) +
    xlim(c(min_x,max_x)) +
    theme_genes()  + theme(axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank()) +
    labs(title=title) + theme(plot.title = element_text(size=16, face="bold",vjust=0.5)) +
    theme(legend.position = "top",legend.text=element_text(size=10),legend.key.size=unit(0.2,"inches"))+guides(color = guide_legend(ncol = NULL, byrow = TRUE) + theme(legend.spacing.x=unit(3,'cm'),legend.border.width=20))

}
## plot gRNA
gscore = read.table(gscore_file,header=T,sep='\t',check.names = F,quote = "")
gscore$direction = ifelse(gscore$StrandOn == "+", 1, 0)
on_score_list = sort(gscore$RuleSet3Score, decreasing = T)
off_sscore_list = sort(gscore$max_CFD_score, decreasing = F)
gscore_label_exon = subset(gscore, grepl("exon", Region))
gscore <- subset(gscore,select = -c(Region,transcript,index))
gscore =  gscore  %>% distinct()
#print (gscore)

### 画所有的gRNA
#print (gscore)
p2 = ggplot(gscore, aes(xmin = StartPos, xmax = StartPos+23, y= gene, fill=StrandOn,forward = direction)) +
	geom_gene_arrow() + labs(x=NULL,y=NULL) + xlim(c(min_x,max_x)) + theme_genes()  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank()) + theme(legend.position = "none")

### 画gRNA对应的位置
exhibit_num = round(dim(gscore)[1]*0.2+0.5) ## 获取整个gRNA得分在前20%的gRNA
#gscore_label = gscore[gscore$RuleSet3Score >= on_score_list[exhibit_num] & gscore$max_CFD_score <= off_sscore_list[exhibit_num],]
gscore_label_exon <- subset(gscore_label_exon,select = -c(Region,transcript,index))
gscore_label_exon = gscore_label_exon %>% distinct()
gscore_label_exon = gscore_label_exon[order(-gscore_label_exon$RuleSet3Score,gscore_label_exon$max_CFD_score),]
gscore_label_exon$label = paste(gscore_label_exon$gRNA_id,gscore_label_exon$gc,gscore_label_exon$RuleSet3Score,gscore_label_exon$mismatch_num,gscore_label_exon$max_CFD_score,sep=':')
if (dim(gscore_label_exon)[1] > 20){ ## 选取前十展示
	gscore_label_exon = gscore_label_exon[1:20,]
}
gscore_label_exon = gscore_label_exon %>% 
	mutate (color = case_when(
		RuleSet3Score < 0.3 ~ "low", ## low
		max_CFD_score >=0.4 ~ "low",
		((RuleSet3Score >= 0.3) & (RuleSet3Score < 0.7) & (max_CFD_score<0.4)) ~ "medium", ## medium
		((RuleSet3Score >= 0.7) & (max_CFD_score<=0.2)) ~ 'high')) ## high

p3_y_min = (dim(gscore_label_exon)[1] - 1) * -0.1 - 0.2
gscore_label_exon['lim_y'] = seq(from = p3_y_min, to=-0.2, by = 0.1)
p3 = ggplot(gscore_label_exon,aes(x=StartPos,y=lim_y,fill=gRNA_id,color=color))+geom_point() +
    theme_genes() +
    ylim(c(p3_y_min-0.2,-0.2)) + theme_classic() + 
	scale_color_manual(values = c("low" = "black", "medium" = "green", "high" = "red")) +
	#scale_color_manual(values = gscore_label_exon$color) +
    theme(legend.position = "bottom",legend.text=element_text(size=10)) +
    theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank()) +
    labs(x=NULL,y=NULL) +
    #geom_text_repel(aes(label=label),size=2,color='black',segment.size=0.2,nudge_x=100,hjust=0)+
	geom_text_repel(aes(label=label),size=2,color='black',direction='y') + 
    xlim(c(min_x,max_x)) + theme(legend.position = "none")

#write.table(gscore_label_exon,file='tmp.txt',sep='\t')
#print (gscore_label_exon)
#ggarrange(p1,p2,p3,ncol=1)
pall = ggarrange(p1,p2,p3,ncol=1) + theme(plot.margin = margin(r=0))
#outfile=paste(outdir,"/",gene,".gRNA_design.png",sep="")
#print (outfile)
#png(outfile,height=600,width=1200)
#print (pall)
#dev.off()
#height = length_transcript
if (length_transcript>30){
    ggsave(file=outfile,plot=pall,height=30,width=10)
}else if (length_transcript>6){
	ggsave(file=outfile,plot=pall,height=ceiling(length_transcript),width=10)
}else{
	ggsave(file=outfile,plot=pall,height=6,width=10)
}
