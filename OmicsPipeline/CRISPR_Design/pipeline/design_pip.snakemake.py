outdir=config['outdir']
indir=config['indir']
script_dir=config['script_dir']
ref=config['ref']
python3=config['python3']
bedtools=config['bedtools']
bowtie2=config['bowtie2']
gpd=config['gpd']
Rscript=config['Rscript']
exon_region=config['exon_region']
max_gc_gRNA=config['max_gc_gRNA']
min_gc_gRNA=config['min_gc_gRNA']
extend=config['extend']
tracrRNA=config['tracrRNA']
primer3_core=config['primer3']

(SAMPLES,) = glob_wildcards("{sample}.txt")

rule all:
    input:
        expand("{outdirs}/{sample}.gRNA.final_primer.csv",sample=SAMPLES,outdirs={outdir})

rule gRNA_design:
    input:
        "{sample}.txt"
    output:
        gRNA_fa = "{outdir}/{sample}.gRNA.fa"
    threads:1
    resources:
        mem_mb = 5000
    shell:
        "{python3} {script_dir}/gRNA_design.py -g {wildcards.sample} -gpd {gpd} --ref {ref} -o {outdir}/{wildcards.sample}.gRNA.on_score.csv -gf {output.gRNA_fa} -a {outdir}/{wildcards.sample}.gene_pred.txt -b {bedtools} --tmp {outdir} --exon {exon_region} -e {extend} -t {tracrRNA} -max_gc {max_gc_gRNA} -min_gc {min_gc_gRNA}"

rule off_target:
    input:
        gRNA_fa = "{outdir}/{sample}.gRNA.fa"
    output:
        "{outdir}/{sample}.gRNA.sam"
    threads:10
    resources:
         mem_mb = 10000
    shell:
        "{bowtie2} -f -N 1 -a -p {threads} --quiet -U {input.gRNA_fa} -x {ref} > {output}"

rule MergeScore:
    input:
        sam = "{outdir}/{sample}.gRNA.sam"
    output:
        mergedscore = "{outdir}/{sample}.on_off.score.csv"
    resources:
        mem_mb = 20000
    shell:
        "{python3} {script_dir}/gRNAmergeScore.py -s {input.sam} -tmp {outdir} -r {ref} -g {wildcards.sample} -o {outdir}/{wildcards.sample}.off_target.csv -m {output.mergedscore} -rg {outdir}/{wildcards.sample}.gRNA.on_score.csv -b {bedtools}"

rule plot:
    input:
        "{outdir}/{sample}.on_off.score.csv"
    output:
        png = "{outdir}/{sample}.gRNA_design.png"
    shell:
        "{Rscript} {script_dir}/plot_Structure_gRNA.R --gene {wildcards.sample} --gpd {outdir}/{wildcards.sample}.gene_pred.txt --gscore {input} --outfile {output.png}"

rule PrimerDesign:
    input:
        "{outdir}/{sample}.on_off.score.csv"
    output:
        primer = "{outdir}/{sample}.gRNA_primer.csv",
        left_primer = "{outdir}/{sample}.primer_left.fa",
        right_primer = "{outdir}/{sample}.primer_right.fa"
    shell:
        "{python3} {script_dir}/gRNA_primer.py --ref {ref} --gRNA {input} --gene {wildcards.sample} --left {output.left_primer} --right {output.right_primer} --outdir {outdir} --bedtools {bedtools} -pr {output.primer} -p {primer3_core}"

rule PrimerSpecific:
    input:
        left_primer = "{outdir}/{sample}.primer_left.fa",
        right_primer = "{outdir}/{sample}.primer_right.fa"
    output:
        "{outdir}/{sample}.primer.sam"
    threads:10
    resources:
        mem_mb = 10000
    shell:
        "{bowtie2} -f --quiet -k 10 -p {threads} --quiet -1 {input.left_primer} -2 {input.right_primer} -x {ref} > {output}"

rule MergePrimer:
    input:
        sam = "{outdir}/{sample}.primer.sam",
        primer = "{outdir}/{sample}.gRNA_primer.csv",
        png = "{outdir}/{sample}.gRNA_design.png"
    output:
        "{outdir}/{sample}.gRNA.final_primer.csv"
    shell:
        "{python3} {script_dir}/gRNAmergePrimer.py --sam {input.sam} --primer {input.primer} --final {output}"
