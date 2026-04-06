#!/usr/bin/env python3
"""
Merges chunked RepeatMasker results back into whole-chromosome coordinates. 
Resolves overlapping regions by finding gaps in annotation coverage within overlap zones and choosing clean cut points.

Called by:
manual execution after 06.2_chunked_repeatmasker.smk
"""


import os
import re
import sys
import argparse
import pandas as pd
from pathlib import Path
import pybedtools
import logging
import shutil

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

def shift_coordinates_gff(gff_file):
    match = re.search(r'chunk_(\d+)\.fasta', gff_file.name)
    offset = int(match.group(1)) if match else 0
    rows = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            rows.append(line.strip().split('\t'))

    if not rows:
        return pd.DataFrame(), None

    df = pd.DataFrame(rows, columns=[
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ])

    true_chromosome_name = df["seqid"][0].rsplit('_', 1)[0]
    df["seqid"] = true_chromosome_name
    df["start"] = df["start"].astype(int) + offset
    df["end"] = df["end"].astype(int) + offset
    df["start"] = df["start"].astype(str)
    df["end"] = df["end"].astype(str)

    return df, true_chromosome_name

def open_rmout(rm_file):
    header = ["SW_score", "perc_div", "perc_del", "perc_ins",
               "query", "q_begin", "q_end", "q_left", "strand",
               "repeat", "repeat_class_family", "r_begin", "r_end", "r_left", "ID", "asterisk"]
    # Read the RepeatMasker output file
    records = []
    with open(rm_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) < 15:
                continue
            rec = dict(zip(header, parts))
            records.append(rec)
    return pd.DataFrame(records)

def shift_coordinates_rm(rm_file):
    match = re.search(r'chunk_(\d+)\.fasta\.out', rm_file.name)
    offset = int(match.group(1)) if match else 0

    rmout = open_rmout(rm_file)

    true_chromosome_name = rmout["query"][0].rsplit('_', 1)[0]
    rmout["query"] = true_chromosome_name
    rmout["q_begin"] = rmout["q_begin"].astype(int) + offset
    rmout["q_end"] = rmout["q_end"].astype(int) + offset
    rmout["q_begin"] = rmout["q_begin"].astype(str)
    rmout["q_end"] = rmout["q_end"].astype(str)
    rmout["q_left"] = "(NA)"

    return rmout, true_chromosome_name



def find_chunk_boundaries(gff_files, chunk_size, true_chromosome_name, chromosome_label, overlap_size, tmp_dir, shifted_dir):
    chunk_starts = [int(re.search(r'chunk_(\d+)', f.name).group(1)) for f in gff_files]
    overlap_starts = chunk_starts[1:]
    overlap_ends = [x + overlap_size for x in overlap_starts]

    overlap_bed = pd.DataFrame({
        "chrom": [true_chromosome_name] * len(overlap_starts),
        "start": overlap_starts,
        "end": overlap_ends
    })
    overlap_bed_path = tmp_dir / f"{chromosome_label}_overlap.bed"
    overlap_bed.to_csv(overlap_bed_path, sep="\t", index=False, header=False)

    bedtools = pybedtools.BedTool(overlap_bed_path)
    shifted_gff_files = [shifted_dir / gff.name for gff in gff_files]
    merged_gff_path = tmp_dir / f"{chromosome_label}_merged_shifted.gff"
    with open(merged_gff_path, "w") as outfile:
        for fname in shifted_gff_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    coverage = bedtools.coverage(str(merged_gff_path), d=True)
    coverage_df = coverage.to_dataframe()
    coverage_df.columns = ["chrom", "start", "end", "overlap_pos", "coverage"]
    coverage_df = coverage_df[coverage_df["coverage"] == 0]
    coverage_df = coverage_df.groupby(["chrom", "start"]).first().reset_index()
    coverage_df["absolute_position"] = coverage_df["start"] + coverage_df["overlap_pos"]

    chunk_info = pd.DataFrame({
        "chrom": [true_chromosome_name] * len(chunk_starts),
        "start": chunk_starts,
        "end": [x + chunk_size for x in chunk_starts]
    })
    chunk_info = chunk_info.merge(coverage_df[["end", "absolute_position"]], on="end", how="left")
    chunk_info["absolute_position"].fillna(chunk_info["end"], inplace=True)
    chunk_info.rename(columns={"absolute_position": "new_end"}, inplace=True)
    chunk_info["new_start"] = chunk_info["new_end"].shift(1)
    chunk_info["new_start"].fillna(0, inplace=True)

    return chunk_info

def add_rm_header(rm_file):
    header = """
        SW   perc perc perc  query                    position in query               matching                        repeat                                                   position in repeat
 score   div. del. ins.  sequence                 begin    end           (left)   repeat                          class/family                                         begin   end    (left)      ID"""

    with open(rm_file, "r") as infile, open(rm_file.with_suffix(".header.out"), "w") as outfile:
        outfile.write(header)
        # add a newline after the header
        outfile.write("\n")
        for line in infile:
            outfile.write(line)

def main(chunk_dir, out_dir=None, overlap_size=100000, chunk_size=100000000, cleanup=True):
    chunk_dir = Path(chunk_dir)
    species = chunk_dir.parts[-2]
    chromosome_label = chunk_dir.parts[-1]
    default_out_dir = Path(f"data/merged_chunks/{species}")
    out_dir = Path(out_dir) if out_dir else default_out_dir
    os.makedirs(out_dir, exist_ok=True)

    gff_files = sorted(chunk_dir.glob("chunk_*.fasta.out.gff"), key=lambda f: int(re.search(r'chunk_(\d+)', f.name).group(1)))
    rm_files = sorted(chunk_dir.glob("chunk_*.fasta.out"), key=lambda f: int(re.search(r'chunk_(\d+)', f.name).group(1)))

    if not gff_files:
        logging.error("No GFF files found.")
        sys.exit(1)

    if not rm_files:
        logging.error("No RepeatMasker .out files found.")
        sys.exit(1)

    logging.info(f"Processing {len(gff_files)} GFF files from {chunk_dir}")
    tmp_dir = out_dir / f"tmp_{chromosome_label}"
    shifted_dir = tmp_dir / "shifted_coordinates"
    shifted_dir.mkdir(parents=True, exist_ok=True)

    for gff in gff_files:
        shifted_gff, true_chromosome_name = shift_coordinates_gff(gff)
        shifted_gff.to_csv(shifted_dir / gff.name, sep="\t", index=False, header=False)

    for rm_file in rm_files:
        shifted_rm, true_chromosome_name = shift_coordinates_rm(rm_file)
        shifted_rm.to_csv(shifted_dir / rm_file.name, sep="\t", index=False, header=False)
    
    chunk_info = find_chunk_boundaries(gff_files, chunk_size, true_chromosome_name, 
                                       chromosome_label, overlap_size, tmp_dir, shifted_dir)

    cropped_dir = tmp_dir / "cropped"
    cropped_dir.mkdir(parents=True, exist_ok=True)

    for gff in gff_files:
        shifted_gff = pd.read_csv(shifted_dir / gff.name, sep="\t", header=None)
        shifted_gff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        chunk_gff_start = int(re.search(r'chunk_(\d+)', gff.name).group(1))
        chunk_row = chunk_info[chunk_info["start"] == chunk_gff_start]
        if chunk_row.empty:
            continue
        new_start = int(chunk_row["new_start"].values[0])
        new_end = int(chunk_row["new_end"].values[0])
        cropped_gff = shifted_gff[
            (shifted_gff["start"].astype(int) >= new_start) &
            (shifted_gff["end"].astype(int) <= new_end)
        ]
        cropped_gff.to_csv(cropped_dir / gff.name, sep="\t", index=False, header=False)

    final_gff_path = out_dir / f"{chromosome_label}.out.gff"
    with open(final_gff_path, "w") as outfile:
        for fname in cropped_dir.glob("*.gff"):
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    for rm_file in rm_files:
        shifted_rm = pd.read_csv(shifted_dir / rm_file.name, sep="\t", header=None)
        shifted_rm.columns = ["SW_score", "perc_div", "perc_del", "perc_ins",
                              "query", "q_begin", "q_end", "q_left", "strand",
                              "repeat", "repeat_class_family", "r_begin", "r_end", "r_left", "ID", "asterisk"]
        chunk_rm_start = int(re.search(r'chunk_(\d+)', rm_file.name).group(1))
        chunk_row = chunk_info[chunk_info["start"] == chunk_rm_start]
        if chunk_row.empty:
            continue
        new_start = int(chunk_row["new_start"].values[0])
        new_end = int(chunk_row["new_end"].values[0])
        cropped_rm = shifted_rm[
            (shifted_rm["q_begin"].astype(int) >= new_start) &
            (shifted_rm["q_end"].astype(int) <= new_end)
        ]
        cropped_rm.to_csv(cropped_dir / rm_file.name, sep="\t", index=False, header=False)

    final_rm_path = Path(f"data/merged_chunks/{species}/{chromosome_label}.out")
    final_rm_path.parent.mkdir(parents=True, exist_ok=True)
    rm_chunks = list(cropped_dir.glob("*.out"))
    rm_chunks.sort(key=lambda x: int(re.search(r'chunk_(\d+)', x.name).group(1))) # Sort by chunk number
    with open(final_rm_path, "w") as outfile: 
        for rm_chunk in rm_chunks: # concatenate each chunk
            with open(rm_chunk) as infile:
                for line in infile:
                    outfile.write(line)

    add_rm_header(final_rm_path)

    logging.info(f"Final merged RepeatMasker .out saved to: {final_rm_path}")

    if cleanup:
        shutil.rmtree(tmp_dir)
        logging.info(f"Temporary files cleaned up from {tmp_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge and adjust GFF chunks with overlaps.")
    parser.add_argument("chunk_dir", help="Directory with chunked GFF files")
    parser.add_argument("--out-dir", default=None, help="Optional output directory")
    parser.add_argument("--overlap-size", type=int, default=100000, help="Overlap size in base pairs")
    parser.add_argument("--chunk-size", type=int, default=100000000, help="Chunk size in base pairs")
    parser.add_argument("--no-cleanup", dest="cleanup", action="store_false", help="Do not clean up intermediate files")
    parser.set_defaults(cleanup=True)
    args = parser.parse_args()
    main(args.chunk_dir, out_dir=args.out_dir, overlap_size=args.overlap_size, chunk_size=args.chunk_size, cleanup=args.cleanup)
