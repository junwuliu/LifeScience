#!/usr/bin/env python3
import pandas as pd
import pickle
import numpy as np
import re
import os
import subprocess
import sys
from mybase.CFDscore import calc_cfd
from mybase.process import processSAM
import argparse

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
import warnings
warnings.filterwarnings('ignore')

parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-s','--sam',help='gRNA alignment sam file',dest='sam',type=str,required=True)
parser.add_argument('-tmp','--tmpdir',help='tmpdir for restore bed file ',dest='tmpdir',type=str,required=False,default='./')
parser.add_argument('-r','--ref',help='genome reference file',dest='ref',type=str,required=True)
parser.add_argument('-g','--gene',help='gene name',dest='gene',type=str,required=True)
parser.add_argument('-b','--bedtools',help='bedtools path',dest='bedtools',type=str,required=False,default='/public/software/bedtools/bedtools')
parser.add_argument('-o','--off',help='all off target site info',dest='off',type=str,required=True)
parser.add_argument('-m','--merge',help='final merged on-score/off-score file',dest='merge',type=str,required=True)
parser.add_argument('-rg','--grna',help='gRNA on score file',dest='grna',type=str,required=True)
args=parser.parse_args()

def get_off_site(align,bedtools,fasta,tmpbed):
	align['POS_START'] = align['POS'] -1
	align['POS_END'] = align['POS_START'] + 23
	align['STRAND'] = align['FLAG'].apply(lambda x:"-" if (int(x) & 0x10) else "+") ## 位运算，匹配0x10为负链
	align_bed = align[['RNAME','POS_START','POS_END','QNAME','SEQ','STRAND']]
	align_bed.to_csv(tmpbed,sep='\t',header=0,index=0)
	bedtools_command = [
		bedtools,
		"getfasta",
		"-fi", fasta,
		"-bed", tmpbed,
		"-name",
		"-s",
	]
	process = subprocess.run(bedtools_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if process.returncode == 0:
	# 命令执行成功，获取stdout中的结果
		result = process.stdout
	#print(result)
	else:
	# 命令执行失败，打印stderr中的错误信息
		print("Command failed with error:")
		print(process.stderr)
		sys.exit(1)
	return result,align_bed

def get_off_info(result):
	multi_result = result.rstrip().split('\n')
	total_off = []
	gene_list = []
	off_target_site = []
	for i in range(0,len(multi_result)-1,2):
		ids = multi_result[i]
		genename,gap,chr,pos = ids.split(":")
	#print (gene)
		start_pos = pos.split("-",maxsplit=1)[0]
		genename = genename.strip('>')
		dna_seq = multi_result[i+1].upper()
		total_off.append({'gRNA_id':genename,'off-target_site':dna_seq,'off_chr':chr,'off_start_pos':start_pos})
	total_off_pd = pd.DataFrame(total_off)
	return total_off_pd

def cal_CFD(seq,off): ## 计算CFDscore
	wt = seq
	pam = off[-2:]
	sg = off[:-3]
	cdf_score = round(calc_cfd(wt,sg,pam),3)
	return cdf_score

def merge(raw_gRNA_file,off_target_alignment,off_file):
	on_score = pd.read_csv(raw_gRNA_file,header=0,sep='\t',dtype=str)
	merge = pd.merge(on_score,off_target_alignment,how='left',left_on='gRNA_id',right_on='gRNA_id')
	multi_alignment = merge[~((merge['chr'] == merge['off_chr']) & (merge['StartPos'] == merge['off_start_pos']))]
	multi_alignment.to_csv(off_file,index=0,sep='\t')
	for index,row in on_score.iterrows():
		gRNA_id = row['gRNA_id']
		gRNA_chr = row['chr']
		gRNA_pos = row['StartPos']
		off_target = off_target_alignment[off_target_alignment['gRNA_id']==gRNA_id]
		if (off_target.shape[0] == 1): ## no off target site
			mismatch_num = 0
			max_CFD_score = 0
			median_CFD_score = 0
		else:
			off_target = off_target[~((off_target['off_chr'] == gRNA_chr) & (off_target['off_start_pos'] == gRNA_pos))]
			mismatch_num = off_target.shape[0]
			max_CFD_score = off_target['CFD_score'].max()
			median_CFD_score = off_target['CFD_score'].median()
		on_score.loc[index,'mismatch_num'] = mismatch_num
		on_score.loc[index,'max_CFD_score'] = max_CFD_score
		on_score.loc[index,'median_CFD_score'] = median_CFD_score
	on_score['mismatch_num'] = on_score['mismatch_num'].astype(int)
	return on_score

if __name__ == '__main__':
	sam = args.sam
	tmpdir = args.tmpdir
	ref = args.ref
	raw_gRNA = args.grna
	bedtools = args.bedtools
	gene = args.gene
	off_file = args.off
	merged_file = args.merge
	tmpbed = args.tmpdir + '/' + gene + '.whole_genome.align.bed'
	whole_alignment = processSAM(sam)
	whole_align_fa,whole_align_bed = get_off_site(whole_alignment,bedtools,ref,tmpbed)
	total_off_pd = get_off_info(whole_align_fa)
	whole_align_bed['POS_START'] = whole_align_bed['POS_START'].astype('str')
	total_off_pd['off_start_pos'] = total_off_pd['off_start_pos'].astype('str')
	off_target_alignment = pd.merge(whole_align_bed,total_off_pd,how='left',left_on=['QNAME','RNAME','POS_START'],right_on=['gRNA_id','off_chr','off_start_pos'])
	off_target_alignment['CFD_score'] = off_target_alignment.apply(lambda x:cal_CFD(x['SEQ'],x['off-target_site']),axis=1)
	merge_score = merge(raw_gRNA,off_target_alignment,off_file)
	merge_score.to_csv(merged_file,sep='\t',index=0)
	os.system('rm -rf {}'.format(tmpbed))
