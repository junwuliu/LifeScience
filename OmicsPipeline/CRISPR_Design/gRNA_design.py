#python3
import argparse
import pandas as pd
import sys
import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import parse
import subprocess
import re
from mybase.calculate_gc import calculate_gc_content
from mybase.process import getfasta
import warnings
warnings.filterwarnings('ignore')

from rs3.seq import predict_seq



parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-g','--gene',help='gene name',dest='gene',type=str,required=True)
parser.add_argument('-tmp','--tmpdir',help='tmpdir for restore bed file ',dest='tmpdir',type=str,required=False,default='./')
parser.add_argument('-gpd','--genepred',help='gene annotation file, GenePred format',dest='gpd',type=str,required=True)
parser.add_argument('-r','--ref',help='genome reference file',dest='ref',type=str,required=True)
parser.add_argument('-gf','--gfasta',help='gRNA fasta file',dest='gfasta',type=str,required=True)
parser.add_argument('-o','--out',help='gene gRNA designed file',dest='outfile',type=str,required=True)
parser.add_argument('-exon','--exon',help='only exon region or not,default is full transcript',dest='exon',type=int,required=False,default=1)
parser.add_argument('-a','--anno',help='gene specifc annotation file, for plot',dest='anno',type=str,required=True)
parser.add_argument('-e','--extend',help='extend transcript uptream and downstream bp for design',dest='extend',type=int,required=False,default=50)
parser.add_argument('-b','--bedtools',help='bedtools executable path',type=str,required=False,default='/public/software/bedtools/bedtools')
parser.add_argument('-t','--tracr',help='tracr RNA name, Hsu2013 or Chen2013',type=str,required=False,default='Hsu2013')
parser.add_argument('-max_gc','--max_gc',help='max gRNA gc theshold',dest='maxgc',required=False,type=int,default=60)
parser.add_argument('-min_gc','--min_gc',help='min gRNA gc theshold',dest='mingc',required=False,type=int,default=40)

args=parser.parse_args()

gene = args.gene
gpd_file = args.gpd
fasta_file = args.ref
outfile = args.outfile
#outexonfile= args.outexonfile
tmpdir = args.tmpdir
gfasta = args.gfasta
tmpbed = tmpdir + "/" + gene + '.bed'
anno = args.anno
max_gc_content = args.maxgc
min_gc_content = args.mingc
extend = args.extend
bedtools = args.bedtools
tracr = args.tracr
## get gene loci
genepred = pd.read_csv(gpd_file,header=None,sep='\t')
genepred.columns = ['transcript_id','chr','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score','geneName','cdsStartStat','cdsEndStat','IstringexonFrames']
gene_annotation = genepred[genepred['geneName']==gene]
gene_annotation.to_csv(anno,index=0,sep='\t')
gene_annotation['gene_transcript'] = gene_annotation['geneName'] + ":" + gene_annotation['transcript_id']
gene_annotation_bed = gene_annotation[['chr','txStart','txEnd','gene_transcript','score','strand']]
gene_annotation_bed['txStartExtend'] = gene_annotation_bed['txStart'].apply(lambda x:x-extend if (x>extend) else x) ## 拓展30bp
gene_annotation_bed['txEndExtend'] = gene_annotation_bed['txEnd'].apply(lambda x:x+extend) ## 如果基因位于染色体尾部？
gene_annotation_bed = gene_annotation_bed[['chr','txStartExtend','txEndExtend','gene_transcript','score','strand']]
strand = gene_annotation_bed.iloc[0,5]
gene_annotation_bed.to_csv(tmpbed,sep='\t',header=0,index=0)

result = getfasta(fasta_file,tmpbed,False,bedtools)

def location(pos,t_gpd):
	location = 'Unknown'
	exon_s = t_gpd.loc[0,'exonStarts'].rstrip(',').split(',')
	#print (exon_s)
	exon_s = list(map(int, exon_s))
	exon_e = t_gpd.loc[0,'exonEnds'].rstrip(',').split(',')
	exon_e = list(map(int, exon_e))
	exon_nums = len(exon_s)
	cds_s = int(t_gpd.loc[0,'cdsStart'])
	cds_e = int(t_gpd.loc[0,'cdsEnd'])
	t_strand = t_gpd.loc[0,'strand']
	if ((pos + 23 <= cds_s) and (pos >= exon_s[0]) ):
		location = "5'UTR" if t_strand == '+' else "3'UTR"
	elif ((pos + 23 > cds_e) and (pos + 23 <= exon_e[-1])):
		location = "3'UTR" if t_strand == '+' else "5'UTR"
	elif (pos < exon_s[0] or pos > exon_e[-1]):
		location = 'intergenic'
	else:
		for m in range(0,exon_nums):
			if (exon_s[m] <= pos + 20 <= exon_e[m]):
				region = 'exon'
				exon_number = m+1 if t_strand == '+' else (exon_nums -m)
				location = 'exon' + "_" + str(exon_number)
			else:
				continue
		if location == 'Unknown':
			location = 'intron'
	return location

### 从序列中寻找NGG特征的gRNA特征
multi_result = result.split('\n')
total_gene = []
for i in range(0,len(multi_result)-1,2):
    id = multi_result[i]
    gene,transcript,gap,chr,pos = id.split(":")
    #print (transcript)
    transcript_gpd = gene_annotation[gene_annotation['transcript_id'] == transcript].reset_index()
    start_pos,end_pos = pos.split("-")
    gene = gene.strip('>')
    forward_seq = multi_result[i+1].upper() ## 获得基因组正链上序列
    dna_seq_obj = Seq(forward_seq)
    reverse_seq = str(dna_seq_obj.reverse_complement()) ## 获得反向互补序列
    matches = re.finditer(r'\w{25}GG\w{3}',forward_seq) ## 寻找正链上30bp长度的gRNA序列(NGG)，用于Rule Set 2 score
    matches_forward = [(match.group(), match.start()) for match in matches] #返回index
    matches = re.finditer(r'\w{25}GG\w{3}',reverse_seq) ## 寻找负链上30bp长度的gRNA序列
    matches_reverse = [(match.group(), match.start()) for match in matches]
    forward_gRNA = pd.DataFrame(matches_forward,columns=['NNNNgRNANNN','index'])
    forward_gRNA['gRNA'] = forward_gRNA['NNNNgRNANNN'].apply(lambda x:x[4:27]) ## 20bp gRNA + PAM
    forward_gRNA['StrandOn'] = "+"
    forward_gRNA['StartPos'] = forward_gRNA['index'].apply(lambda x:int(x)+int(start_pos) + 4) ## gRNA的起始
    reverse_gRNA =pd.DataFrame(matches_reverse,columns=['NNNNgRNANNN','index'])
    reverse_gRNA['gRNA'] = reverse_gRNA['NNNNgRNANNN'].apply(lambda x:x[4:27])
    reverse_gRNA['StrandOn'] = "-"
    reverse_gRNA['StartPos'] = reverse_gRNA['index'].apply(lambda x:int(end_pos)-1-int(x)-30 +4) ## 和取出的长度保持一致
    gene_gRNA = pd.concat([forward_gRNA,reverse_gRNA],axis=0)
    gene_gRNA['Region'] = gene_gRNA['StartPos'].apply(lambda x:location(x,transcript_gpd)) ## 获得gRNA处于intergenic/exon/intron区间
    gene_gRNA.insert(loc=0,column='chr',value=chr)
    gene_gRNA.insert(loc=0,column='transcript',value=transcript)
    gene_gRNA.insert(loc=0,column='gene',value=gene)
    total_gene.append(gene_gRNA)
total_gene_merge = pd.concat(total_gene,axis=0)
total_gene_merge['gc'] = total_gene_merge['gRNA'].apply(lambda x:calculate_gc_content(x))
total_gene_merge = total_gene_merge[(total_gene_merge['gc'] <= max_gc_content) & (total_gene_merge['gc'] >= min_gc_content) ]
#total_gene_merge.to_csv(outfile,sep='\t',index=0)
print (total_gene_merge.shape)
if (args.exon == 1):
	total_gene_merge = total_gene_merge[total_gene_merge['Region'].str.contains('exon')] #### 暂时只找exon
	print (total_gene_merge.shape)
else:
	pass

# 获得全局gRNA fasta file,用于脱靶评估
unique_gRNA = total_gene_merge['gRNA'].unique().tolist()
outgFasta=open(gfasta,'w')
unique_id = []
for i in range(0,len(unique_gRNA)):
    gRNA_id = 'gRNA_' + str(i)
    outgFasta.writelines('>'+gRNA_id +'\n')
    outgFasta.writelines(unique_gRNA[i] +'\n')
    unique_id.append({'gRNA':unique_gRNA[i],'gRNA_id':gRNA_id})

unique_id_pd = pd.DataFrame(unique_id)
total_gene_merge = pd.merge(total_gene_merge,unique_id_pd,how='left',left_on='gRNA',right_on='gRNA')
total_gene_merge['RuleSet3Score'] = pd.Series(predict_seq(total_gene_merge['NNNNgRNANNN'].to_list(), sequence_tracr=tracr))
total_gene_merge['RuleSet3Score'] = total_gene_merge['RuleSet3Score'].round(3)
total_gene_merge.to_csv(outfile,sep='\t',index=0)

#gene_exon_gRNA = total_gene_merge[total_gene_merge['Region'].str.contains('exon')]
#gene_exon_gRNA.to_csv(outexonfile,sep='\t',index=0)
os.system('rm -rf {}'.format(tmpbed))

