import pandas as pd
import subprocess
import sys
import argparse
from mybase.process import getfasta

parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-b','--bedtools',help='bedtools path',dest='bedtools',type=str,required=False,default='/public/software/bedtools/bedtools')
parser.add_argument('-p','--primer3',help='primer3 path',dest='primer3',type=str,required=False,default='/public/software/primer3-main/src/primer3_core')
parser.add_argument('-o','--outdir',help='output dir',dest='outdir',type=str,required=True)
parser.add_argument('-r','--ref',help='reference fasta path',dest='reference',type=str,required=True)
parser.add_argument('-g','--gRNA',help='gRNA score file',dest='grna',type=str,required=True)
parser.add_argument('-G','--gene',help='gene name',dest='gene',type=str,required=True)
parser.add_argument('-pr','--presult',help='primer result',dest='presult',type=str,required=True)
parser.add_argument('-le','--left',help='left primer fasta file',dest='left',type=str,required=True)
parser.add_argument('-ri','--right',help='right primer fasta file',dest='right',type=str,required=True)
args=parser.parse_args()

primer3_templete = '''PRIMER_PRODUCT_SIZE_RANGE=150-300
PRIMER_EXPLAIN_FLAG=1
PRIMER_NUM_RETURN=5
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=17
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60
PRIMER_MIN_TM=58
PRIMER_MAX_TM=62
PRIMER_MIN_GC=45
PRIMER_MAX_GC=55
PRIMER_MAX_POLY_X=4
='''

def write_templete(templete,templete_seq,grna_id):
	outfile = open(templete,'w')
	outfile.writelines("SEQUENCE_TEMPLATE="+templete_seq+'\n')
	outfile.writelines("SEQUENCE_ID="+grna_id+'\n')
	outfile.writelines(primer3_templete)
	outfile.close()

def run_primer3(gene_primer_templete):
	primer3_command = [
		primer3_core, 
		gene_primer_templete,
	]
	process = subprocess.run(primer3_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if process.returncode == 0:
	# 命令执行成功，获取stdout中的结果
		result = process.stdout
		#print(result)
	else:
		# 命令执行失败，打印stderr中的错误信息
		print("Command failed with error:")
		print(process.stderr)
		result = "NotFound"
		sys.exit(1)
	return result
	


if __name__ == '__main__':
	bedtools = args.bedtools
	primer3_core = args.primer3
	reference = args.reference
	grna_score_file = args.grna
	outdir = args.outdir
	gene = args.gene
	total_primer_result = args.presult
	left_primer_fa_file = args.left
	right_primer_fa_file = args.right
	
	gRNA_score = pd.read_csv(grna_score_file,header=0,sep='\t')
	gRNA_score = gRNA_score[['gene','chr','gRNA','StrandOn','StartPos','gc','gRNA_id','RuleSet3Score','mismatch_num','max_CFD_score','median_CFD_score']].drop_duplicates()
	gRNA_score = gRNA_score.drop_duplicates() ## 去重
	gRNA_score['TempleteStart'] = gRNA_score['StartPos']  - 150 +23  ## 前后拓展150bp
	gRNA_score['TempleteEnd'] = gRNA_score['StartPos'] + 150
	templete_bed = gRNA_score[['chr','TempleteStart','TempleteEnd','gRNA_id','gene','StrandOn']]
	tmpbed = outdir + '/' + gene + ".primer_templete.bed"
	templete_bed.to_csv(tmpbed,sep='\t',header=0,index=0)
	templete_fasta = getfasta(reference,tmpbed,False,bedtools)
	multi_result = templete_fasta.rstrip().split('\n')
	## 运行primer3 和 处理primer3 结果文件
	total_primer_info = []
	for i in range(0,len(multi_result)-1,2):
		ids = multi_result[i]
		grna_id,gap,chr,pos = ids.split(":")
		grna_id = grna_id.replace(">","")
		templete_seq = multi_result[i+1]
		gene_primer_templete = "".join([outdir,"/",gene,".primer.config"])
		write_templete(gene_primer_templete,templete_seq,grna_id)
		primer_result = run_primer3(gene_primer_templete)
		primer_results = primer_result.split('\n')
		primer_results = primer_results[25:-2] ## 该数值取决于templte多长 
		num=0
		for j in range(0,len(primer_results)-1,23): #依次获取primer信息
			num = num +1
			primer_info = primer_results[j:j+23]
			primer_info = pd.DataFrame(primer_info,columns=['info'])
			primer_info[['flag','name']] = primer_info['info'].str.split("=",expand=True)
			primer_info = primer_info.drop('info',axis=1)
			primer_info['flag'] = primer_info['flag'].str.replace("_{}".format(num-1),'')
			primer_info['gRNA_id'] = grna_id
			primer_info = primer_info.pivot(index='gRNA_id',columns='flag', values='name').reset_index()
			primer_info.insert(loc=1,column='PrimerNum',value=num)
			total_primer_info.append(primer_info)
	total_primer = pd.concat(total_primer_info,axis=0)
	total_primer = pd.concat(total_primer_info,axis=0)
	total_primer= total_primer[['gRNA_id','PrimerNum','PRIMER_PAIR_PRODUCT_SIZE',
	  'PRIMER_LEFT_SEQUENCE','PRIMER_RIGHT_SEQUENCE',
	  'PRIMER_LEFT_GC_PERCENT','PRIMER_RIGHT_GC_PERCENT',
	  'PRIMER_LEFT_TM','PRIMER_RIGHT_TM','PRIMER_LEFT','PRIMER_RIGHT',
	  'PRIMER_PAIR_PRODUCT_TM','PRIMER_PAIR_PENALTY',
	  ]]
	total_primer['PRIMER_LEFT_SEQUENCE'] = total_primer['PRIMER_LEFT_SEQUENCE'].str.upper()
	total_primer['PRIMER_RIGHT_SEQUENCE'] = total_primer['PRIMER_RIGHT_SEQUENCE'].str.upper()
	total_primer[['Primer_Left_start','Primer_Left_length']] = total_primer.PRIMER_LEFT.str.split(',',expand=True)
	total_primer[['Primer_Right_start','Primer_Right_length']] = total_primer.PRIMER_RIGHT.str.split(',',expand=True)
	## 处理primer3 结果文件完成
	gRNA_score_combine = gRNA_score[['chr','StartPos','TempleteStart','TempleteEnd','gRNA_id','gene','StrandOn']].drop_duplicates().reset_index(drop=True)
	total_primer = pd.merge(total_primer,gRNA_score_combine,how='left',left_on='gRNA_id',right_on='gRNA_id')
	total_primer['left_primer_start'] = total_primer['TempleteStart'].astype('int') + total_primer['Primer_Left_start'].astype('int')
	total_primer['right_primer_end'] = total_primer['TempleteStart'].astype('int') + total_primer['Primer_Right_start'].astype('int') + 1
	total_primer.to_csv(total_primer_result,sep='\t',index=0)
	## 生成primer的对称fasta文件
	left_primer_fa = open(left_primer_fa_file,'w')
	right_primer_fa = open(right_primer_fa_file,'w')
	for index,row in total_primer.iterrows():
		left_id = ":".join([row['gRNA_id'],str(row['PrimerNum'])])
		right_id = ":".join([row['gRNA_id'],str(row['PrimerNum'])])				   
		left_primer_fa.writelines(">"+left_id+'\n'+row['PRIMER_LEFT_SEQUENCE']+'\n')
		right_primer_fa.writelines(">"+right_id+'\n'+row['PRIMER_RIGHT_SEQUENCE']+'\n')
	left_primer_fa.close()
	right_primer_fa.close()

