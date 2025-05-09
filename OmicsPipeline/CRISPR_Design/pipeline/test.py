import os
#configfile: "design_pip.yaml"
indir=config['indir']
outdir=config['outdir']

#exec("print(indir)")
exec("os.chdir(indir)")
#exec("print(os.getcwd())")

#(SAMPLES,) = glob_wildcards(f"{config['indir']} + "*.txt"")
#SAMPLES = glob_wildcards("config['indir']/*.txt")
(SAMPLES,) = glob_wildcards("{sample}.txt")
exec("print (SAMPLES)")
#exec("print (AAA)")
#input_files = expand("{indir}/{sample}.txt", sample=SAMPLES,indir=indir)
#output_files = expand("{indir}/{sample}.csv", sample=SAMPLES,indir=indir)
rule all:
	input:
		expand("{indir}/{sample}.csv", sample=SAMPLES,indir=indir)

rule bwa_map:
	input:
		"{sample}.txt"
	output:
		"{sample}.csv"
	shell:
		"echo {input} > {output}"
