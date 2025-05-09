(SAMPLES,)= glob_wildcards("{sample}.txt")
rule all:  
    input:
        expand('{sample_name}.1_1_fastqc.zip',sample_name=SAMPLES)
rule fastqc:
    input:
        fq='{sample_name}.1_1.fastq.gz'
    output:
        '{sample_name}.1_1_fastqc.zip'
    log:
        '{sample_name}.1_1.log'
    params:
        outdir='qc'
    shell:
        'fastqc {input[fq]} -o {params[outdir]} 1>{log[0]} 2>&1'
