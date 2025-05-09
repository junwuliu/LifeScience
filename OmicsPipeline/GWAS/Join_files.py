#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
import argparse
import pandas as pd
import sys
import numpy as np
import os

parser=argparse.ArgumentParser(description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i1','--input1',help='input1 file',dest='input1',required=True)
parser.add_argument('-i2','--input2',help='input2 file',dest='input2',required=True)
#parser.add_argument('-c1','--input1_char',help='input1 join char',dest='char1',required=True)
#parser.add_argument('-c2','--input2_char',help='input2 join char',dest='char2',required=True)
parser.add_argument('-o','--output',help='output file',dest='output',required=True)
args=parser.parse_args()

data1 = pd.read_csv(args.input1,header=0,sep=r'\s+')
data2 = pd.read_csv(args.input2,header=0,sep=r'\s+')

#print (data1.head())
#print (data2.head())

joint_data = pd.merge(data1,data2,how='left',left_on=["CHR","SNP"],right_on=["CHR","SNP"])
print (joint_data.head())

joint_data.to_csv(args.output,index=0,sep='\t')
