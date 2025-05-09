#! /usr/bin/env python3
# -*- coding: UTF-8 -*-
import argparse
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

parser=argparse.ArgumentParser(description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i','--input',help='input file',dest='input',required=True)
parser.add_argument('-o','--output',help='output file',dest='output',required=True)
args=parser.parse_args()

data = pd.read_csv(args.input,header=None,sep=' ')

data = data[[1,2,3]]
data.columns = ['IID','PCA1','PCA2']

plt.scatter(data['PCA1'], data['PCA2'], c='green')
plt.legend(title='labels',bbox_to_anchor=(1.005, 0.5),loc=3,borderaxespad=0)

plt.savefig(args.output)
plt.show()
