import os
import pandas as pd
import glob

excluded_species = {}
chromosomes_list = []

with open('data/chromosome_list_provisional.txt') as f:
    for line in f:
        key = line.strip()
        if not key:
            continue
        sample_name = key.split('/')[0]
        if sample_name in excluded_species:
            continue
        chromosomes_list.append(key)

rule all:
    input:
        expand("data/repeatmasker_filtered/{chromosome}.fna.out", chromosome=chromosomes_list),
        expand("data/solo_ratio/{chromosome}.solo_list", chromosome=chromosomes_list)


rule run_repeatmasker:
    input:
        reference = "data/chromosomes/{chromosome}.fna"
    params:
        # base dir of chromosome
        sample = lambda w: w.chromosome.split('/')[0]
    output:
        out = "data/repeatmasker/{chromosome}.fna.out",
        out_gff = "data/repeatmasker/{chromosome}.fna.out.gff",
    resources:
        memgb="64gb",
        walltime="4:00:00",
        cpus=8
    #group: "10min"
    shell:
        """
        module add repeatmasker
        export TMPDIR=$SCRATCHDIR
        export OMP_NUM_THREADS=$PBS_NUM_PPN
        mkdir -p data/repeatmasker/{params.sample}
        RepeatMasker -pa $PBS_NUM_PPN -lib data/ltr_library/{params.sample}/TE_DLplus.fasta -gff -dir data/repeatmasker/{params.sample} {input.reference} -no_is -nolow
        """

rule filter_repeatmasker:
    input:
        rmout = "data/repeatmasker/{chromosome}.fna.out"
    output:
        filtered = "data/repeatmasker_filtered/{chromosome}.fna.out"
    params:
        # base dir of chromosome
        sample = lambda w: w.chromosome.split('/')[0]
    resources:
        memgb = "8gb",
        walltime = "1:00:00",
        cpus = 1
    group: "1min"
    shell:
        """
        module add mambaforge-22.9.0
        conda activate /storage/plzen1/home/kratka/.conda/envs/merge_chunks
        mkdir -p data/repeatmasker_filtered/{params.sample}
        python workflow/scripts/filter_repeatmasker_records.py {input.rmout} > {output.filtered}
        """

rule solo_finder:
    input:
        "data/repeatmasker_filtered/{chromosome}.fna.out"
    params:
        # base dir of chromosome
        sample = lambda w: w.chromosome.split('/')[0]
    output:
        "data/solo_ratio/{chromosome}.solo_list"
    resources:
        memgb="16gb",
        walltime="4:00:00",
        cpus=1
    group: "10s"
    shell:
        """  
        export OMP_NUM_THREADS=$PBS_NUM_PPN
        export TMPDIR=$SCRATCHDIR
        perl workflow/scripts/solo_finder.pl -i {input} -info data/ltr_library/{params.sample}/lib.LTR.info > {output}
        awk '{{print $1,$2,$3}}' {output} > data/solo_ratio/{wildcards.chromosome}.solo.bed
        """

