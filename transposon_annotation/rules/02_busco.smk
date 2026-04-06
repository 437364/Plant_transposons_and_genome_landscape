import os
import pandas as pd

# create a dictionary from input table of samples that contains name + sra id
samples_information = pd.read_csv("data/accessions_genome.tsv", sep='\t', index_col=False)
sample_id = list(samples_information['Sample_id'])
ftp_link = list(samples_information['FTP_link'])
samples_dict = dict(zip(sample_id, ftp_link))

rule all:
    input:
        expand("data/genomes/{sample}.fna.gz", sample=samples_dict.keys()),
        expand("data/chromosomes/{sample}_ok", sample=samples_dict.keys()),
        expand("data/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt", 
               sample=samples_dict.keys(), lineage=config.get("busco_lineage", "viridiplantae_odb10"))
        #expand("data/assembly_stats/{sample}_gaps.tsv", sample=samples_dict.keys()),
        #"data/assembly_stats/all_samples_stats.tsv"

rule run_busco:
    input:
        fna = "data/genomes/{sample}.fna.gz"
    output:
        summary = "data/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt",
        full = "data/busco/{sample}/run_{lineage}/full_table.tsv"
    params:
        lineage = config.get("busco_lineage", "viridiplantae_odb10"),  # default plant lineage
        outdir = "data/busco/{sample}",
        mode = "genome",
    resources:
        memgb = "128gb",
        walltime = "24:00:00",
        cpus = 8
    shell:
        """
        set +u
        # Uncompress genome if needed
        mkdir -p {params.outdir}
        gunzip -c {input.fna} > ${{SCRATCH}}/{wildcards.sample}_temp.fna
        
        # Run BUSCO
        module add mambaforge
        conda activate busco
        busco \
            -i ${{SCRATCH}}/{wildcards.sample}_temp.fna \
            -o {wildcards.sample} \
            -l {params.lineage} \
            -m {params.mode} \
            -c {resources.cpus} \
            --out_path data/busco/ -f
        
        # Clean up temp file
        rm ${{SCRATCH}}/{wildcards.sample}_temp.fna
        set -u
        """