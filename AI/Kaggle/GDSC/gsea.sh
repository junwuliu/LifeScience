#! /bin/sh
expression=$1
outdir=$2

#python3 /public/home/liujunwu/workdir/DataSets/Kaggle/GDSC/transfer_cellLine_ssGSEA.py -i $1 -o $1.tmp -e /public/home/liujunwu/workdir/DataSets/Kaggle/GDSC/GSEA/header
/public/home/liujunwu/software/miniconda3/envs/R4/bin/gseapy ssgsea -d $1.tmp -g /public/home/liujunwu/workdir/scripts/GNN_Reactome/GTEx_exp/GSEA/ReactomePathways.GSEA.gmt -o $2 --min-size 0 --max-size 5000 --sample-norm rank -c zscore -p 10
