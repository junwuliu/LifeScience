ref=GRCh38.p14_genomic.fna
gpd=NCBI.GRCh38.annotation.gpd
#bedtools=`which bedtools`
#primer3=`which primer3_core`
#Rscript=`which Rscript`
bedtools=/public/software/bedtools/bedtools
primer3=/public/software/primer3-main/src/primer3_core
Rscript=/public/software/R/bin/Rscript
tmpdir=./

python3 gRNA_design.py -g PDCD1 -gpd $gpd --ref $ref -o PDCD1.gRNA.on_score.csv -gf PDCD1.gRNA.fa -a PDCD1.genepred.txt -b $bedtools --exon 1 --tmpdir $tmpdir

bowtie2 -f -N 1 -a -p 10 --quiet -U PDCD1.gRNA.fa -x $ref > PDCD1.gRNA.sam

python3 gRNAmergeScore.py -s PDCD1.gRNA.sam -tmp ./ -r $ref -g PDCD1 -o PDCD1.off_target.csv -m PDCD1.on_off.score.csv -rg PDCD1.gRNA.on_score.csv -b $bedtools

$Rscript plot_Structure_gRNA.R --gene PDCD1 --gpd PDCD1.genepred.txt --gscore PDCD1.on_off.score.csv --outfile PDCD1.gRNA_design.png

python3 gRNA_primer.py --ref $ref --gRNA PDCD1.on_off.score.csv --gene PDCD1 --left PDCD1.primer_left.fa --right PDCD1.primer_right.fa --outdir ./ --bedtools $bedtools -pr PDCD1.gRNA_primer.csv -p $primer3

bowtie2 -f --quiet -k 10 -p 10 --quiet -1 PDCD1.primer_left.fa -2 PDCD1.primer_right.fa -x $ref > PDCD1.primer.sam

python3 gRNAmergePrimer.py --sam PDCD1.primer.sam --primer PDCD1.gRNA_primer.csv --final PDCD1.gRNA.final_primer.csv
