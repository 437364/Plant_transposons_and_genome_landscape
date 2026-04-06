import glob

chunk_files = glob.glob("data/chromosomes_chunked/*/*/chunk_*.fasta")
chunks = [
    f.replace("data/chromosomes_chunked/", "").replace(".fasta", "")
    for f in chunk_files
]

rule all:
    input:
        expand("data/repeatmasker_chunked/{chunk}.fasta.out", chunk=chunks)

rule run_repeatmasker_chunk:
    input:
        fasta = "data/chromosomes_chunked/{chunk}.fasta"
    params:
        sample = lambda w: w.chunk.split('/')[0]
    output:
        out = "data/repeatmasker_chunked/{chunk}.fasta.out",
        gff = "data/repeatmasker_chunked/{chunk}.fasta.out.gff"
    resources:
        memgb = "128gb",
        walltime = "168:00:00",
        cpus = 16
    shell:
        """
        module add repeatmasker
        export TMPDIR=$SCRATCHDIR
        export OMP_NUM_THREADS=$PBS_NUM_PPN
        mkdir -p $(dirname {output.out})
        RepeatMasker -pa {resources.cpus} \
                     -lib data/ltr_library/{params.sample}/TE_DLplus.fasta \
                     -gff \
                     -dir $(dirname {output.out}) \
                     {input.fasta} \
                     -no_is -nolow
        """