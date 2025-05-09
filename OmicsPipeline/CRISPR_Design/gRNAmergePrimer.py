import pandas as pd
import sys
import argparse

from mybase.process import processSAM

parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-s','--sam',help='primer alignment sam file',dest='sam',type=str,required=True)
parser.add_argument('-p','--primer',help='gRNA primer result',dest='primer',type=str,required=True)
parser.add_argument('-f','--final',help='gRNA final primer result, merged specific result',dest='final',type=str,required=True)
args=parser.parse_args()


total_primer = pd.read_csv(args.primer,header=0,sep='\t',dtype=str)
whole_alignment = processSAM(args.sam)
whole_alignment = whole_alignment[whole_alignment['RNEXT'] == '=']
primer_count = whole_alignment['QNAME'].value_counts().reset_index()
primer_count['multi_align'] = primer_count['count'].apply(lambda x:'Yes' if x>2 else "No")
primer_count.columns = ['JointName','count','multi_align']

total_primer['JointName'] = total_primer['gRNA_id'] + ":" + total_primer['PrimerNum'].astype(str)
total_primer = pd.merge(total_primer,primer_count,how='left',left_on='JointName',right_on='JointName')
total_primer['PosRelatedPrimer'] = total_primer.apply(lambda x:'inside' if ((int(x['StartPos']) > int(x['left_primer_start'])+int(x['Primer_Left_length'])) & (int(x['StartPos']) +23 < int(x['right_primer_end']) -int(x['Primer_Right_length']) )) else 'outside' if ((int(x['StartPos']) < int(x['left_primer_start'])) & (int(x['StartPos']) >  int(x['right_primer_end']))) else 'overlap',axis=1)
total_primer.to_csv(args.final,index=0,sep='\t')
