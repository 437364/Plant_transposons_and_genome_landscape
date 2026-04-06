import os
import pandas as pd
import glob

chromosomes_list = []
species_with_large_chromosomes = ["Hordeum_vulgare", "Secale_cereale", "Triticum_aestivum", "Triticum_monococcum"]

for filename in glob.glob('data/chromosomes/**/*.fna.gz', recursive=True):
    # extract "sample_name/chr_name" as key
    key = filename.split('/')[2] + '/' + filename.split('/')[3].split('.')[0]
    species = key.split('/')[0]
    if species in species_with_large_chromosomes:
        chromosomes_list.append(key)

CHUNK_SIZE = 100_000_000
OVERLAP = 100_000
rule all:
    input:
        expand("data/chromosomes_chunked/{chromosome}.is.chunked", chromosome=chromosomes_list)

rule chunk_chromosome:
    input:
        fasta = "data/chromosomes/{chromosome}.fna.gz"
    output:
        done = "data/chromosomes_chunked/{chromosome}.is.chunked"
    params:
        outdir = "data/chromosomes_chunked/{chromosome}",
        chunk_size = CHUNK_SIZE,
        overlap = OVERLAP
    resources:
        memgb="16gb",
        walltime="24:00:00",
        cpus=1
    group: "1min"
    shell:
        """
        module add python36-modules
        mkdir -p {params.outdir}
        python workflow/scripts/chunk_fasta.py {input.fasta} {params.outdir} {params.chunk_size} {params.overlap}
        touch {output.done}
        """