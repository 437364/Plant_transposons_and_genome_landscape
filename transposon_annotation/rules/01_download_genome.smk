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
        expand("data/reports/{sample}_assembly_report.txt", sample=samples_dict.keys())
        #expand("data/chromosomes/{sample}_ok", sample=samples_dict.keys())
    
rule download_genome:
    params:
        ftp_link = lambda w: samples_dict[w.sample],
        accession = lambda w: samples_dict[w.sample].split("/")[-1]
    output:
        fna = "data/genomes/{sample}.fna.gz"
    resources:
        memgb="8gb",
        walltime="24:00:00",
        cpus=1
    group: "1min"
    run:
        print(wildcards.sample)
        # ftp link
        print(f'{params.ftp_link}/{params.accession}_genomic.fna.gz')
        
        shell("""
        #export OMP_NUM_THREADS=$PBS_NUM_PPN
        mkdir -p data/genomes
        wget {params.ftp_link}/{params.accession}_genomic.fna.gz -O {output.fna}
        """)

rule download_assembly_report:
    params:
        ftp_link = lambda w: samples_dict[w.sample],
        accession = lambda w: samples_dict[w.sample].split("/")[-1]
    output:
        report = "data/reports/{sample}_assembly_report.txt"
    resources:
        memgb="8gb",
        walltime="24:00:00",
        cpus=1
    group: "1min"
    shell:
        """
        mkdir -p data/reports
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 {params.ftp_link}/{params.accession}_assembly_report.txt -O {output.report}
        """


rule extract_chromosomes:
    input:
        fna = "data/genomes/{sample}.fna.gz",
        report = "data/reports/{sample}_assembly_report.txt"
    output:
        flag = "data/chromosomes/{sample}_ok"
    resources:
        memgb="8gb",
        walltime="24:00:00",
        cpus=1
    group: "5min"
    run:
        import pandas as pd
        from pathlib import Path
        import subprocess

        report = pd.read_csv(input.report, sep='\t', comment='#', header=None)
        report.columns = [
            "Sequence-Name",
            "Sequence-Role",
            "Assigned-Molecule",
            "Assigned-Molecule-Location/Type",
            "GenBank-Accn",
            "Relationship",
            "RefSeq-Accn",
            "Assembly-Unit",
            "Sequence-Length",
            "UCSC-style-name"
        ]

        # Keep only chromosomes
        chrom_df = report[
            (report["Assigned-Molecule-Location/Type"] == "Chromosome") &
            (report["Sequence-Role"] == "assembled-molecule")
        ].copy()

        # Use RefSeq-Accn if available, else fall back to GenBank-Accn
        is_missing = chrom_df["RefSeq-Accn"].str.lower().eq("na") | chrom_df["RefSeq-Accn"].isna()
        chrom_df["Extraction-ID"] = chrom_df["RefSeq-Accn"].where(~is_missing, chrom_df["GenBank-Accn"])

        # Prepend "chr" if sequence name is only digits
        acc_to_name = {}
        for acc, name in zip(chrom_df["Extraction-ID"], chrom_df["Sequence-Name"]):
            if str(name).isdigit():
                name = f"chr{name}"
            acc_to_name[acc] = name

        # Index fasta if needed
        fna_path = Path(input.fna).with_suffix('')  # uncompressed path

        if not fna_path.exists():
            subprocess.run(f"gunzip -c {input.fna} > {fna_path}", shell=True, check=True)

        if not fna_path.with_suffix('.fai').exists():
            subprocess.run(f"module add samtools && samtools faidx {fna_path}", shell=True, check=True)

        outdir = Path(f"data/chromosomes/{wildcards.sample}")
        outdir.mkdir(parents=True, exist_ok=True)

        for acc, name in acc_to_name.items():
            outpath = outdir / f"{name}.fna"
            # Extract the sequence, keep original header
            cmd = f"module add samtools && samtools faidx {fna_path} {acc} > {outpath}"
            subprocess.run(cmd, shell=True, check=True)

        Path(output.flag).touch()
