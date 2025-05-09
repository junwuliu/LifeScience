目录文件说明：
├── gRNA_design.py                 gRNA设计主程序
├── gRNAmergeScore.py              合并gRNA on-score和off-score的程序
├── gRNA_primer.py                 gRNA引物设计主程序, primer3 模板相关参数设定在该代码内部
├── gRNAmergePrimer.py             合并gRNA引物结果的程序
├── mybase
│   ├── calculate_gc.py				gRNA的GC含量计算子程序
│   ├── CFDscore.py					CFD计算代码
│   ├── mismatch_score.pkl			CFD得分评估模型
│   ├── pam_scores.pkl				CFD PAM区得分评估模型
│   ├── process.py					SAM文件处理/bedtools取系列子程序
├── pipeline
│   ├── design_pip.snakemake.py		snakemake流程代码
│   ├── design_pip.yaml				snakemake配置文件
│   └── run.sh						运行snakemake的代码
├── plot_Structure_gRNA.R           画基因结构和gRNA位置得分的代码
├──	NCBI.GRCh38.annotation.gpd		NCBI hg38基因注释文件
├── 逻辑图.png						gRNA设计和primer设计逻辑图
└── run.sh							单样本运行范例


参考基因组为NCBI hg38(版本：GCF_000001405.40_GRCh38.p14)，由于文件太大, 需自行前往NCBI官网下载并建立索引。

结果为三部分：
1、gRNA设计结果：gene.on_off.score.csv,记录了gRNA对应的转录本信息，位置信息，gc含量，on-score和off-score的得分
2、gRNA位置图：gRNA_design.png，详细含义已向徐博说明。
3、gRNA primer设计结果：gene.gRNA.final_primer.csv, 记录了gRNA最优的5条引物信息，multi_align代表是否在基因组其他位置有比对(特异性问题)，PosRelatedPrimer代表gRNA和primer的相对位置。详细含义已向徐博说明。
其余文件均为中间文件，结果生成后可以删除。

运行环境：
python3 3.9
bowtie2 2.2.5
R	4.1.2
Primer3	2.6.1
bedtools v2.30.0
samtools 1.14

python3依赖包：
pandas	2.1.1
argparse 1.1
numpy	1.23.1
Bio	1.79
re	2.2.1
rs3	0.0.16

R模块：
ggpubr
ggrepel
patchwork
ggplot2
gggenes
gridExtra
tidyr
dplyr



