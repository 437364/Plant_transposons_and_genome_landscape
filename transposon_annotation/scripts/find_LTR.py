#!/usr/bin/env python3
"""
Replacement for LTR_retriever's find_LTR.pl. 
Extracts LTR boundary coordinates from the DANTE-LTR library by matching TE sequences back to GFF3 annotations. 
Faster than BLAST-based approach.

Called by:
05_prepare_library.smk
"""

import argparse
from Bio import SeqIO
import pandas as pd

def extract_te_and_ltr_annotations(annotations):
    te_annotations = annotations[annotations['feature_type'] == 'transposable_element'].copy()
    te_annotations["ID"] = te_annotations["attributes"].str.extract(r'ID=(\w+)')

    ltr_annotations = annotations[annotations['feature_type'] == 'long_terminal_repeat'].copy()
    ltr_annotations["Parent"] = ltr_annotations.attributes.str.extract(r'Parent=(\w+)')
    ltr_annotations["Side"] = ltr_annotations.attributes.str.extract(r'LTR=(\w+)')
    return te_annotations, ltr_annotations

def process_sequence(record, te_annotations, ltr_annotations):
    chr, start, end = record.id.split("#")[0].rsplit("_",2)
    start, end = int(start), int(end)

    element_row = te_annotations[(te_annotations["chr"] == chr) & 
                                 (te_annotations["start"] == start) & 
                                 (te_annotations["end"] == end)]
    
    if element_row.empty:
        return
    
    element_id = element_row["ID"].values[0]
    ltrs = ltr_annotations[ltr_annotations["Parent"] == element_id]

    if len(ltrs) != 2:
        return
    
    ltr5 = ltrs[ltrs["Side"] == "5LTR"]
    ltr3 = ltrs[ltrs["Side"] == "3LTR"]
    
    if ltr5.empty or ltr3.empty:
        return
    
    ltr5_len = ltr5["end"].values[0] - ltr5["start"].values[0] + 1
    ltr3_len = ltr3["end"].values[0] - ltr3["start"].values[0] + 1

    ltr3_start = len(record.seq) - ltr3_len + 1
    ltr5_info = f"{record.id}\t1\t{ltr5_len}\t{ltr5_len}"
    ltr3_info = f"{record.id}\t{ltr3_start}\t{len(record.seq)}\t{ltr3_len}"

    print(ltr5_info)
    print(ltr3_info)

def main(lib, gff):
    gff_columns = ['chr', 'source', 'feature_type', 'start', 'end', 'strand', 'score', 'phase', 'attributes']
    annotations = pd.read_csv(gff, sep='\t', header=None, names=gff_columns, comment='#')
    te_annotations, ltr_annotations = extract_te_and_ltr_annotations(annotations)
    library = SeqIO.to_dict(SeqIO.parse(lib, "fasta"))
    for record in library.values():
        process_sequence(record, te_annotations, ltr_annotations)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        This script prints lib.LTR.info to the standard output, 
        describing the LTRs found in the input library.  
    """)
    parser.add_argument("lib", help="Input FASTA library containing full-length elements.")
    parser.add_argument("gff", help="Input GFF file with full-length element annotations.")
    
    args = parser.parse_args()
    main(args.lib, args.gff)