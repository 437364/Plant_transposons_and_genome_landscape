import os
import pandas as pd
import glob

chromosomes_list = []
for filename in glob.glob('data/chromosomes/**/*.fna', recursive=True):
    # extract "sample_name/chr_name" as key
    relative_path = filename.replace('data/chromosomes/', '').replace('.fna', '')
    chromosomes_list.append(relative_path)


rule all:
    input:
        expand("data/ltr/{chromosome}.tes.gff3", chromosome=chromosomes_list)


rule run_dante:
    input:
        fna = "data/chromosomes/{chromosome}.fna"
    output:
        domains = "data/ltr/{chromosome}.domains.gff3"
    resources:
        memgb="32gb",
        walltime="24:00:00",
        cpus=4
    group: "5min"
    shell:
        """
        module add mambaforge-22.9.0
        conda activate /storage/plzen1/home/kratka/.conda/envs/dante_ltr
        export TMPDIR=$SCRATCHDIR
        export OMP_NUM_THREADS=$PBS_NUM_PPN

        dante -q {input.fna} -o {output.domains} -c $PBS_NUM_PPN
        """
    
rule run_dante_ltr:
    input:
        unzip = "data/chromosomes/{chromosome}.fna",
        domains = "data/ltr/{chromosome}.domains.gff3"
    output:
        ltr = "data/ltr/{chromosome}.tes.gff3"
    resources:
        memgb="32gb",
        walltime="48:00:00",
        cpus=1
    group: "10min"
    shell:
        """
        module add mambaforge-22.9.0
        conda activate /storage/plzen1/home/kratka/.conda/envs/dante_ltr
        export TMPDIR=$SCRATCHDIR
        export OMP_NUM_THREADS=$PBS_NUM_PPN
        dante_ltr -g {input.domains} -s {input.unzip} -o data/ltr/{wildcards.chromosome}.tes -M 3 -c $PBS_NUM_PPN
        """
