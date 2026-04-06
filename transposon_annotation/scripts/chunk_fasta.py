#!/usr/bin/env python
"""
Script to chunk a FASTA file into smaller pieces with optional overlap.
Usage:
python chunk_fasta.py input.fasta.gz output_dir chunk_size overlap
- input.fasta.gz: Input FASTA file (can be gzipped)
- output_dir: Directory to save the chunked FASTA files
- chunk_size: Size of each chunk in bases
- overlap: Number of bases to overlap between chunks
Example:
python chunk_fasta.py genome.fasta.gz chunks/ 1000000 10000


Called by:
06.1_chunk_large_chromosomes.smk
"""


import sys
import gzip
from pathlib import Path
from Bio import SeqIO

fasta_path = sys.argv[1]
outdir = Path(sys.argv[2])
chunk_size = int(sys.argv[3])
overlap = int(sys.argv[4])

outdir.mkdir(parents=True, exist_ok=True)

# Read gzipped or plain FASTA
open_fn = gzip.open if fasta_path.endswith(".gz") else open
with open_fn(fasta_path, "rt") as f:
    for record in SeqIO.parse(f, "fasta"):
        seq = record.seq
        length = len(seq)
        for start in range(0, length, chunk_size - overlap):
            end = min(start + chunk_size, length)
            chunk_seq = seq[start:end]
            chunk_name = f"chunk_{start}.fasta"
            with open(outdir / chunk_name, "w") as out_f:
                out_f.write(f">{record.id}_{start}-{end}\n{chunk_seq}\n")