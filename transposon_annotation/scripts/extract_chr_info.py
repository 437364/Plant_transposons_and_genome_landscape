#!/usr/bin/env python3
"""
Used to build:
chromosome_sizes.tsv

Extract chromosome information from FASTA files and format as TSV rows.
Processes a directory of .fna files and outputs species, chromosome ID, and length.

"""

import sys
from pathlib import Path


def get_fasta_length(fasta_file):
    """Calculate total sequence length from a FASTA file."""
    length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines and empty lines
            if line and not line.startswith('>'):
                length += len(line)
    return length


def extract_chr_id_from_filename(filename):
    """Extract chromosome ID from filename (removes .fna extension)."""
    return filename.replace('.fna', '')


def process_species_directory(species_dir, species_name):
    """Process all .fna files in a species directory and yield TSV rows."""
    species_path = Path(species_dir)
    
    if not species_path.exists():
        print(f"Warning: Directory {species_dir} does not exist", file=sys.stderr)
        return
    
    # Find all .fna files and sort them
    fasta_files = sorted(species_path.glob('*.fna'))
    
    if not fasta_files:
        print(f"Warning: No .fna files found in {species_dir}", file=sys.stderr)
        return
    
    for fasta_file in fasta_files:
        chr_id = extract_chr_id_from_filename(fasta_file.name)
        length = get_fasta_length(fasta_file)
        yield f"{species_name}\t{chr_id}\t{length}"


def main():
    if len(sys.argv) < 3:
        print("Usage: python extract_chr_info.py <species_name> <directory_path> [<species_name> <directory_path> ...]")
        print("\nExample: python extract_chr_info.py Cuscuta_australis data/chromosomes/Cuscuta_australis")
        sys.exit(1)
    
    # Process pairs of (species_name, directory_path)
    args = sys.argv[1:]
    for i in range(0, len(args), 2):
        if i + 1 >= len(args):
            print(f"Warning: Missing directory path for {args[i]}", file=sys.stderr)
            break
        
        species_name = args[i]
        directory_path = args[i + 1]
        
        for row in process_species_directory(directory_path, species_name):
            print(row)


if __name__ == '__main__':
    main()