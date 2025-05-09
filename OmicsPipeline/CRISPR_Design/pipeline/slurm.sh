#! /bin/sh
#snakemake -np --configfile design_pip.yaml -s design_pip.snakemake.py
#snakemake --dag -s ../design_pip.snakemake.py --configfile ../design_pip.yaml --forceall |dot -Tpdf > dag.pdf
snakemake --cluster "sbatch --mem={resources.mem_mb} --cpus-per-task={threads}" --jobs 10 --configfile design_pip.yaml -s design_pip.snakemake.py --cores 10
