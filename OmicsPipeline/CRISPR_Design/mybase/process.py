import pandas as pd
import subprocess
import sys


def processSAM(sam_file):
    alignment = pd.read_csv(sam_file,header=None,sep='\t',comment="@",names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','AS','XS','XN','XM','XO','XG','NM','MD','YS','YT:Z'])
    #alignment = alignment.iloc[:,0:11]
    alignment = alignment[['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL']]
    #alignment.columns = ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','AS','XS','XN','XM','XO','XG','NM','MD','YT:Z']
    #alignment.columns = ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL']
    alignment = alignment[~alignment['RNAME'].str.contains('NT_|NW_')]
    return alignment

def getfasta(ref,bed,strand,bedtools):
	bedtools_command = [
		bedtools,
		"getfasta", 
		"-fi", ref,
		"-bed", bed,
		"-name",
]
	if strand: ## 如果需要根据正负链取序列
		bedtools_command.append("-s")
	else:
		pass
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
	return result

