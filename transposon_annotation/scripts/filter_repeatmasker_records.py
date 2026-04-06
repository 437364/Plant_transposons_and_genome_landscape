#!/usr/bin/env python3
"""
Filter RepeatMasker .out rows with inconsistent library sequence sizes.
Based on this issue with RepeatMasker output:
https://github.com/Dfam-consortium/RepeatMasker/issues/357

Called by:
06_find_elements.smk

Usage:
python filter_repeatmasker_records.py input.out -o output.tsv --criterion both --keep-header
- input.out: RepeatMasker .out file
- output.tsv: Output TSV file (optional, default: stdout)
- --criterion: Filtering rule (default: both)
    'end' = drop rows where repeat_end > library_size
    'sum' = drop rows where (repeat_end + repeat_left|begin) != library_size
    'both' (default) = apply both filters
- --keep-header: Prepend a simple header line with column names to the TSV
Example:
python filter_repeatmasker_records.py genome.out -o filtered.tsv --criterion end --keep-header


"""


import argparse
import sys
import re
import pandas as pd

COLUMNS = [
    "SW_score","perc_div","perc_del","perc_ins","query_seq",
    "query_start","query_end","query_left","strand",
    "matching_repeat","repeat_class","repeat_begin","repeat_end","repeat_left","ID"
]

def parse_args():
    p = argparse.ArgumentParser(
        description="Filter RepeatMasker .out rows with inconsistent library sequence sizes."
    )
    p.add_argument("rmout", help="RepeatMasker .out file")
    p.add_argument("-o", "--output", help="Output TSV file (default: stdout)")
    p.add_argument(
        "--criterion",
        choices=["end", "sum", "both"],
        default="both",
        help=("Filtering rule: "
              "'end' = drop rows where repeat_end > library_size; "
              "'sum' = drop rows where (repeat_end + repeat_left|begin) != library_size; "
              "'both' (default) = apply both filters.")
    )
    p.add_argument(
        "--keep-header",
        action="store_true",
        help="Prepend a simple header line with column names to the TSV."
    )
    return p.parse_args()

def process_repeatmasker_output(rmout: str) -> pd.DataFrame:
    # Read after the 3-line header
    repeatmasker_data = []
    with open(rmout, "r") as f:
        lines = []
        for line in f:
            if re.match(r"^\s*\d", line):
                lines.append(line)

    for line in lines:
        line = line.strip()
        if not line:
            continue
        parts = re.split(r"\s+", line, maxsplit=14)
        if len(parts) >= 15:
            repeatmasker_data.append(parts[:15])

    df = pd.DataFrame(repeatmasker_data, columns=COLUMNS)

    # Make numeric copies for calculations; keep originals intact for output
    # Handle parentheses in repeat_begin/repeat_end/repeat_left
    def to_int_series(s):
        return s.str.replace(r"[()]", "", regex=True).astype(int)

    df["_repeat_begin_i"] = to_int_series(df["repeat_begin"])
    df["_repeat_end_i"]   = to_int_series(df["repeat_end"])
    df["_repeat_left_i"]  = to_int_series(df["repeat_left"])

    # Library size from matching_repeat "X_Y_Z" -> Z - Y + 1
    def lib_size_from_matching(m):
        try:
            parts = m.split("_")
            end = int(parts[-1]); start = int(parts[-2])
            return end - start + 1
        except Exception:
            return pd.NA

    df["_lib_size"] = df["matching_repeat"].apply(lib_size_from_matching).astype("Int64")

    # rm_size inferred from RM columns depends on strand
    # + strand: repeat_end + repeat_left == lib_size
    # C strand: repeat_end + repeat_begin == lib_size
    df["_rm_size"] = df["_repeat_end_i"]  # base
    is_plus = df["strand"] == "+"
    is_comp = df["strand"] == "C"

    df.loc[is_plus, "_rm_size"] = df.loc[is_plus, "_repeat_end_i"] + df.loc[is_plus, "_repeat_left_i"]
    df.loc[is_comp, "_rm_size"] = df.loc[is_comp, "_repeat_end_i"] + df.loc[is_comp, "_repeat_begin_i"]

    # Checks
    df["_end_ok"] = df["_repeat_end_i"] <= df["_lib_size"]
    df["_sum_ok"] = df["_rm_size"] == df["_lib_size"]

    return df

def filter_df(df: pd.DataFrame, criterion: str) -> pd.DataFrame:
    if criterion == "end":
        mask = df["_end_ok"].fillna(False)
    elif criterion == "sum":
        mask = df["_sum_ok"].fillna(False)
    else:  # both
        mask = (df["_end_ok"] & df["_sum_ok"]).fillna(False)

    filtered = df.loc[mask, COLUMNS].copy()
    return filtered, mask

def main():
    args = parse_args()
    df = process_repeatmasker_output(args.rmout)
    filtered, mask = filter_df(df, args.criterion)

    # Diagnostics to stderr
    total = len(df)
    kept = mask.sum()
    dropped = total - kept
    end_fail = (~df["_end_ok"].fillna(True)).sum()
    sum_fail = (~df["_sum_ok"].fillna(True)).sum()

    print(f"Total records: {total}", file=sys.stderr)
    print(f"Kept records:  {kept}", file=sys.stderr)
    print(f"Dropped:       {dropped}", file=sys.stderr)
    print(f"End>size fails: {end_fail}", file=sys.stderr)
    print(f"Sum!=size fails: {sum_fail}", file=sys.stderr)

    # Output TSV
    if args.output:
        with open(args.output, "w") as out:
            if args.keep_header:
                out.write("\t".join(COLUMNS) + "\n")
            filtered.to_csv(out, sep="\t", header=False, index=False)
    else:
        if args.keep_header:
            sys.stdout.write("\t".join(COLUMNS) + "\n")
        filtered.to_csv(sys.stdout, sep="\t", header=False, index=False)

if __name__ == "__main__":
    main()
