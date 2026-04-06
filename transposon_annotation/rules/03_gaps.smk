import os
import pandas as pd
import glob


chromosomes_list = []
for filename in glob.glob('data/chromosomes/**/*.fna', recursive=True):
    # extract "sample_name/chr_name" as key
    relative_path = filename.replace('data/chromosomes/', '').replace('.fna', '')
    chromosomes_list.append(relative_path)

def get_chromosome_file(wildcards):
    """Return the chromosome file, checking for .gz or uncompressed version"""
    base_path = f"data/chromosomes/{wildcards.chromosome}.fna"
    gz_path = f"{base_path}.gz"
    
    if os.path.exists(gz_path):
        return gz_path
    elif os.path.exists(base_path):
        return base_path
    else:
        # Return gz by default, let Snakemake handle the error
        return gz_path

rule all:
    input:
        expand("data/gaps/{chromosome}.gaps.bed", chromosome=chromosomes_list)

rule find_chromosome_gaps:
    input:
        fna = get_chromosome_file
    output:
        gaps_bed = "data/gaps/{chromosome}.gaps.bed"
    params:
        sample = lambda w: w.chromosome.split('/')[0]
    resources:
        memgb = "8gb",
        walltime = "1:00:00",
        cpus = 1
    group: "1min"
    shell:
        """
        module add seqtk
        
        mkdir -p data/gaps/{params.sample}
        seqtk gap {input.fna} > {output.gaps_bed}

        """
