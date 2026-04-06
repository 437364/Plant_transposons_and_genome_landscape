import os
import glob
import pandas as pd

""" # check whether the ltr dir contains all chromosomes
chromosomes_ready = {}
for filename in glob.glob('data/ltr/**/*.tes.gff3', recursive=True):
    # extract "sample_name" "chr_name"
    sample_name = filename.split('/')[2]
    #chromosome_name = filename.split('/')[3].split('.')[0]
    if sample_name not in chromosomes_ready:
        chromosomes_ready[sample_name] = 1
    else:
        chromosomes_ready[sample_name] += 1

samples_information = pd.read_csv("data/accessions_genome.tsv", sep='\t', index_col=False)

for sample in list(samples_information['Sample_id']):
    chromosomes_total = samples_information[samples_information['Sample_id'] == sample]['Chromosome_number'].values[0]
    if sample not in chromosomes_ready:
        print(f"Sample {sample} is missing")
    elif chromosomes_ready[sample] < chromosomes_total:
        print(f"Sample {sample} is missing {chromosomes_total - chromosomes_ready[sample]} chromosomes")
    else:
        print(f"Sample {sample} is ready")
        # create {sample}_ready file
        with open(f"data/ltr/{sample}_ready", 'w') as f:
            f.write("ready")
"""
samples_list = [] 

for filename in glob.glob('data/ltr/*_ready'):
    # extract "sample_name" as key
    key = filename.split('/')[-1].split('_')[0:-1]
    key = '_'.join(key)
    samples_list.append(key)

rule all:
    input:
        #expand("data/ltr_samplewise/{sample}.tes.gff3", sample=samples_list),
        expand("data/ltr_library/{sample}/TE_DLplus.fasta", sample=samples_list),
        expand("data/ltr_library/{sample}/lib.LTR.info", sample=samples_list)
        

rule merge_chromosomes:
    input: 
        "data/ltr/{sample}_ready"
    output:
        merged_gff = "data/ltr_samplewise/{sample}.tes.gff3"
    resources:
        memgb="32gb",
        walltime="4:00:00",
        cpus=1
    group: "1min"
    shell:
        """
        module add mambaforge-22.9.0
        conda activate /storage/plzen1/home/kratka/.conda/envs/dante_ltr
        export TMPDIR=$SCRATCHDIR
        cat data/ltr/{wildcards.sample}/*.tes.gff3 > {output.merged_gff}
        dante_ltr_summary -g {output.merged_gff} -o data/ltr_samplewise/{wildcards.sample}_summary -a
        """

rule prepare_library:
    input:
        "data/ltr_samplewise/{sample}.tes.gff3"
    output:
        library = "data/ltr_library/{sample}/TE_DLplus.fasta",
        reference = temp("data/tmp/{sample}.fna")
    resources:
        memgb="32gb",
        walltime="24:00:00",
        cpus=1
    params:
        library_dir = "data/ltr_library/{sample}"
    group: "5min"
    shell:
        """
        module add mambaforge-22.9.0
        conda activate /storage/plzen1/home/kratka/.conda/envs/dante_ltr
        export TMPDIR=$SCRATCHDIR
        gunzip -c data/genomes/{wildcards.sample}.fna.gz > {output.reference}
        dante_ltr_to_library -g {input} -s {output.reference} -o {params.library_dir} -c $PBS_NUM_PPN 
        """

rule find_LTR:
    input:
        dante_lib = "data/ltr_library/{sample}/TE_DLplus.fasta",
        gff = "data/ltr_samplewise/{sample}.tes.gff3"
    output:
        ltr_info = "data/ltr_library/{sample}/lib.LTR.info"
    resources:
        memgb="32gb",
        walltime="24:00:00",
        cpus=1
    group: "5min"
    shell:
        """
        module add mambaforge-22.9.0
        conda activate /storage/brno2/home/kratka/.conda/envs/python_tools
        export OMP_NUM_THREADS=$PBS_NUM_PPN
        export TMPDIR=$SCRATCHDIR
        python workflow/scripts/find_LTR.py {input.dante_lib} {input.gff} > {output.ltr_info}
        """