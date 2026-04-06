#!/usr/bin/env python3
"""
TE_pipeline.py

Jednoduchý end-to-end pipeline pro jednu TE rodinu:
- parsuje CD-HIT .clstr soubory pro zvolenou family
- vytváří tabulku clusterů (identity, cluster_id, size, centroid, members)
- vytváří growth_df (max_size, avg_size, n_clusters)
- extrahuje centroidy z FASTA (pokud FASTA chybí, automaticky sestaví z jednotlivých cNN.fa)
- spustí MAFFT -> alignment centroids
- spustí FastTree -> strom (pokud dostupný)
- spočte pairwise p-distance z alignu, udělá hierarchické clusterování (sublines)
- vykreslí heatmapu, trajectory plots
- detekuje HGT kandidáty (odlehlé centroidy)
Outputs jsou v --outdir.

Requires: python3, mafft (optional), FastTree (optional)
Python packages: biopython, pandas, numpy, scipy, matplotlib, seaborn

Usage example:
python3 TE_pipeline.py \
  --family Athila \
  --clstr_dir ./Jeff_cd_hit/ \
  --fasta ./Jeff_RTs/Jeff_Athila_RTs.fa \
  --outdir ./Athila_results \
  --threads 8
"""

import argparse, os, re, glob, subprocess, sys, json
from collections import defaultdict, Counter
from math import isnan

# third-party
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from io import StringIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import seaborn as sns

# seaborn theme & palette
heatmap_palette = sns.color_palette("cubehelix_r", as_cmap=True)
sns.set_theme(context="paper", style="white")

# -----------------------
# Helpers
# -----------------------
def run_cmd(cmd, stdout=None, stderr=None):
    print("CMD:", " ".join(cmd))
    res = subprocess.run(cmd, stdout=stdout, stderr=stderr)
    if res.returncode != 0:
        print("Warning: command failed:", " ".join(cmd), file=sys.stderr)
    return res.returncode

def find_clstr_files(clstr_dir, family, spec=None):
    # pattern matches like Jeff_Athila_RTs_c080.fa.clstr  or anyprefix_Athila_RTs_c*.clstr
    pat = os.path.join(clstr_dir, f"*_{family}_RTs_c*.fa.clstr")
    files = sorted(glob.glob(pat))
    if len(files) == 0:
        # sometimes clstr filenames are without .fa before .clstr
        pat2 = os.path.join(clstr_dir, f"*_{family}_RTs_c*.clstr")
        files = sorted(glob.glob(pat2))
    return files

def parse_clstr_file(path):
    import re
    import pandas as pd

    seq_re = re.compile(r">(.*?)...\s")          # ID = vše mezi ">" a 1. mezerou
    id_re = re.compile(r"at\s+([\d\.]+)%")    # identity
    centroid_re = re.compile(r"\*$")          # centroid indicator

    rows = []
    cluster_id = None

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
                continue

            # parse sequence id
            m = seq_re.search(line + " ")     # přidáme 1 mezeru pro jistotu
            if not m:
                continue

            seq_id = m.group(1)

            # remove trailing ... if present
            seq_id = seq_id.replace("...", "")

            # identity
            m2 = id_re.search(line)
            identity = float(m2.group(1)) if m2 else None

            # centroid?
            is_centroid = bool(centroid_re.search(line))

            rows.append({
                "cluster_id": cluster_id,
                "seq_id": seq_id,
                "is_centroid": is_centroid,
                "identity": identity
            })

    df = pd.DataFrame(rows)

    # assign centroid to whole cluster
    df["centroid"] = df.groupby("cluster_id")["seq_id"].transform(
        lambda x: x[df.loc[x.index, "is_centroid"]].iloc[0]
    )

    return df

def build_cluster_table(clstr_files):
    rows = []
    for path in clstr_files:
        fname = os.path.basename(path)
        m = re.search(r"c(\d{2,3})", fname)
        identity = m.group(1) if m else "NA"

        clusters = parse_clstr_file(path)

        # group by cluster_id to get members
        for cluster_id, group in clusters.groupby("cluster_id"):
            members = group["seq_id"].tolist()
            centroid = group["centroid"].iloc[0] if len(group) > 0 else ""
            rows.append({
                "file": fname,
                "identity": identity,
                "cluster_id": f"CL{cluster_id}",
                "size": len(members),
                "centroid": centroid,
                "members": ",".join(members)
            })

    df = pd.DataFrame(rows, columns=["file","identity","cluster_id","size","centroid","members"])
    
    # normalize identity to numeric (e.g., '090' -> 90)
    df["identity_num"] = df["identity"].astype(str).str.extract(r"(\d+)").astype(float)
    df = df.sort_values(["identity_num","cluster_id"], ascending=[True, True])
    
    return df


def ensure_family_fasta(fasta_path, clstr_files, family, outdir, spec):
    """
    If fasta_path provided and exists -> use it.
    Else, try to assemble one from source FASTA names that share family in clstr files.
    """
    if fasta_path and os.path.exists(fasta_path):
        return fasta_path
    # try to find fasta files that match pattern in same directory as clstr files
    # pattern example: Jeff_Athila_RTs.fa
    cldir = os.path.dirname(clstr_files[0])
    # search one level up as well (user structure)
    candidates = []
    for d in [cldir, os.path.join(cldir, ".."), outdir]:
        pat = os.path.join(d, f"*{spec}_{family}*RTs*.fa")
        cand = glob.glob(pat)
        candidates.extend(cand)
    candidates = sorted(set(candidates))
    if len(candidates)>0:
        # if there's single candidate, return it
        print("Found candidate FASTA(s):", candidates)
        # if multiple, try exact match with spec if provided
        return candidates[0]
    # else, try to assemble centroids from small per-identity fasta files (if present)
    # look for any *.fa files sharing _cNN pattern
    assembled = os.path.join(outdir, f"{spec}_{family}_all_from_cN.fa")
    seqs_written = 0
    for cl in clstr_files:
        base = os.path.basename(cl).replace(".clstr","")
        fa_candidate = os.path.join(os.path.dirname(cl), base)
        if os.path.exists(fa_candidate):
            # read fasta and write into assembled
            for rec in SeqIO.parse(fa_candidate, "fasta"):
                SeqIO.write(rec, open(assembled, "a"), "fasta")
                seqs_written += 1
    if seqs_written>0:
        print(f"Assembled FASTA {assembled} from per-identity FASTA files ({seqs_written} records).")
        return assembled
    raise FileNotFoundError("No FASTA provided and could not find or assemble one for family "+family)

def extract_centroids_table(df_clusters, identity_top=None):
    """
    Choose centroids to represent lineages. Strategy:
     - prefer centroids from highest identity (largest identity_num)
     - unique them (one centroid per cluster across levels)
    """
    # order by descending identity
    if identity_top:
        df_sel = df_clusters[df_clusters["identity"]==identity_top]
    else:
        df_sel = df_clusters.sort_values("identity_num", ascending=False)
    # pick unique centroids in order
    seen=set()
    rows=[]
    for _,r in df_sel.iterrows():
        c = r["centroid"]
        if c not in seen and c != "":
            seen.add(c)
            rows.append({"centroid":c, "identity":r["identity"], "size":r["size"]})
    return pd.DataFrame(rows)

def write_fasta_for_ids(ids, src_fasta, out_fa):
    seqs = {}
    for rec in SeqIO.parse(src_fasta, "fasta"):
        seqs[rec.id] = rec
    with open(out_fa, "w") as out:
        n=0
        for id0 in ids:
            # sometimes id in clstr is 'TE85|Jeff|Tork' while fasta id may be 'TE85|Jeff|Tork' or 'TE85'
            if id0 in seqs:
                SeqIO.write(seqs[id0], out, "fasta"); n+=1
            else:
                # try prefix match (id until first '|')
                key = id0.split("|")[0]
                if key in seqs:
                    SeqIO.write(seqs[key], out, "fasta"); n+=1
                else:
                    # try find any record whose id startswith key
                    found=False
                    for k in seqs.keys():
                        if k.startswith(key):
                            SeqIO.write(seqs[k], out, "fasta"); n+=1; found=True; break
                    if not found:
                        print("Warning: centroid", id0, "not found in", src_fasta, file=sys.stderr)
        print(f"Wrote {n} centroid sequences to {out_fa}")
    return out_fa

def compute_pairwise_pdist(aln_file):
    # read alignment via Bio.AlignIO as SeqRecord objects
    aln = AlignIO.read(aln_file, "fasta")
    n = len(aln)
    names = [rec.id for rec in aln]
    L = aln.get_alignment_length()
    mat = np.zeros((n,n), dtype=float)
    for i in range(n):
        si = str(aln[i].seq)
        for j in range(i+1, n):
            sj = str(aln[j].seq)
            # count positions where both not gap
            pos = 0
            diff = 0
            for a,b in zip(si,sj):
                if a == "-" or b == "-":
                    continue
                pos += 1
                if a != b:
                    diff += 1
            pd = diff / pos if pos>0 else 1.0
            mat[i,j] = mat[j,i] = pd
    return names, mat

# -----------------------
# Main
# -----------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--family", required=True, help="TE family name, e.g. Athila")
    p.add_argument("--clstr_dir", required=True, help="directory with .clstr files")
    p.add_argument("--spec", required=True, help="plant species abbreviation")
    p.add_argument("--fasta", default=None, help="family fasta (optional). If not provided pipeline will try to find or assemble one")
    p.add_argument("--outdir", default="./TE_pipeline_out", help="output directory")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--k", type=int, default=5, help="number of sublineages to cut from centroid tree (default 5)")
    args = p.parse_args()

    family = args.family
    clstr_dir = args.clstr_dir
    fasta_path = args.fasta
    outdir = args.outdir
    threads = args.threads
    spec = args.spec
    k = args.k

    os.makedirs(outdir, exist_ok=True)
    # 1) find clstr files
    clstr_files = find_clstr_files(clstr_dir, family)
    
    if len(clstr_files)==0:
        print("No clstr files found for family", family, "in", clstr_dir)
        sys.exit(1)
    print("Found clstr files:", clstr_files)

    # 2) parse cluster tables
    df_clusters = build_cluster_table(clstr_files)
    df_clusters.to_csv(os.path.join(outdir, f"{spec}_{family}_clusters.tsv"), sep="\t", index=False)
    print("Wrote cluster table:", os.path.join(outdir, f"{spec}_{family}_clusters.tsv"))

    # 3) compute growth_df (per identity)
    growth = df_clusters.groupby("identity").agg(
        max_size=("size","max"),
        avg_size=("size","mean"),
        n_clusters=("cluster_id","count")
    ).reset_index()
    # normalize identity order
    growth["identity_num"] = growth["identity"].astype(float)
    growth = growth.sort_values("identity_num")
    growth.to_csv(os.path.join(outdir, f"{spec}_{family}_growth_df.tsv"), sep="\t", index=False)
    print("Wrote growth df:", os.path.join(outdir, f"{spec}_{family}_growth_df.tsv"))

    # 4) plots: trajectory of max_size and avg_size
    plt.figure(figsize=(8,4))
    sns.lineplot(data=growth, x="identity_num", y="max_size", marker="o", label="max_size")
    sns.lineplot(data=growth, x="identity_num", y="avg_size", marker="o", label="avg_size")
    plt.gca().invert_xaxis()
    plt.xlabel("identity (cXX)")
    plt.ylabel("cluster size")
    plt.title(f"{spec} {family} cluster growth")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{spec}_{family}_growth_plot.png"))
    plt.close()
    print("Saved growth plot")

    # heatmap (max_size)
    mat = growth.set_index("identity")["max_size"].to_frame().T
    plt.figure(figsize=(max(4, len(mat.columns)*0.4), 2))
    #sns.heatmap(mat, annot=True, fmt="g")
    sns.heatmap(mat, annot=True, fmt="g", cmap=heatmap_palette)
    
    plt.title(f"{spec} {family} max_size heatmap")
    plt.savefig(os.path.join(outdir, f"{spec}_{family}_maxsize_heatmap.png"))
    plt.close()

    # 5) prepare centroid fasta
    centroids_df = extract_centroids_table(df_clusters)
    centroid_ids = list(centroids_df["centroid"].unique())
    # determine fasta source
    fasta_real = ensure_family_fasta(fasta_path, clstr_files, family, outdir, spec)
    cent_fasta = os.path.join(outdir, f"{spec}_{family}_centroids.fa")
    write_fasta_for_ids(centroid_ids, fasta_real, cent_fasta)

    # 6) MSA via MAFFT
    aln_file = os.path.join(outdir, f"{spec}_{family}_centroids.aln.fa")
    mafft_bin = shutil_which("mafft")
    if mafft_bin:
        cmd = [mafft_bin, "--auto", "--reorder", "--thread", str(threads), cent_fasta]
        with open(aln_file, "w") as outfh:
            run_cmd(cmd, stdout=outfh)
    else:
        print("MAFFT not found in PATH; writing centroid fasta only. Please run MAFFT manually to produce alignment", aln_file)
        # write unaligned copy to aln_file for compatibility
        SeqIO.write(list(SeqIO.parse(cent_fasta, "fasta")), aln_file, "fasta")

    # 7) compute pairwise distances from alignment
    names, distmat = compute_pairwise_pdist(aln_file)
    # save distance matrix
    pd.DataFrame(distmat, index=names, columns=names).to_csv(os.path.join(outdir, f"{spec}_{family}_centroid_pdist.tsv"), sep="\t")
    
    # Load distance matrix as DataFrame if needed (you already have distmat and names)
    df = pd.DataFrame(distmat, index=names, columns=names)

    from matplotlib import gridspec
    from matplotlib.colors import Normalize
    from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

    # seaborn theme
    sns.set_theme(context="notebook", style="white")

    # -----------------------
    # compute linkage (phylogram)
    # -----------------------
    condensed_dist = squareform(distmat)
    Z = linkage(condensed_dist, method="average")

    # order matrix by tree leaves
    order = leaves_list(Z)
    df_ord = df.iloc[order, order]

    # -----------------------
    # dynamic color scaling
    # -----------------------
    vmin = distmat.min()
    vmax = distmat.max()
    margin = 0.05 * (vmax - vmin)
    norm = Normalize(
        vmin=max(0, vmin - margin),
        vmax=min(1, vmax + margin)
    )

    # -----------------------
    # figure layout (NO label column)
    # -----------------------
    fig = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(
        nrows=1,
        ncols=3,
        width_ratios=[3, 14, 1],  # dendrogram : heatmap : colorbar
        wspace=0.05
    )

    ax_dendro = fig.add_subplot(gs[0])
    ax_heat   = fig.add_subplot(gs[1])
    ax_cbar   = fig.add_subplot(gs[2])

    # -----------------------
    # dendrogram (left)
    # -----------------------
    dendrogram(
        Z,
        orientation="left",
        no_labels=True,
        ax=ax_dendro,
        color_threshold=None,
    )
    ax_dendro.invert_yaxis()
    ax_dendro.set_axis_off()

    # -----------------------
    # heatmap (center)
    # -----------------------
    im = ax_heat.imshow(
        df_ord.values,
        cmap=heatmap_palette,
        norm=norm,
        interpolation="nearest",
        aspect="auto",
    )

    # COMPLETELY remove ticks, numbers, frame
    ax_heat.set_axis_off()

    # -----------------------
    # colorbar (right)
    # -----------------------
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.ax.tick_params(labelsize=24, width=2, length=8)
    cbar.set_label(
        "p-distance",
        fontsize=26,
        rotation=270,
        labelpad=25
    )

    # -----------------------
    # title
    # -----------------------
    ax_heat.set_title(
        f"Pairwise distance clustering for {spec} {family} TEs",
        fontsize=18,
        pad=16
    )

    # -----------------------
    # save figure
    # -----------------------
    output_png = os.path.join(outdir, f"{spec}_{family}_centroid_pdist.png")
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Clustered heatmap saved to {output_png}")

    # 8) hierarchical clustering -> sublineages
    # convert to condensed distance for linkage
    if distmat.shape[0] > 1:
        condensed = squareform(distmat)
        Z = linkage(condensed, method="average")
        labels = fcluster(Z, t=k, criterion="maxclust")
    else:
        labels = np.array([1])
    sub_df = pd.DataFrame({"centroid": names, "lineage": labels})
    sub_df.to_csv(os.path.join(outdir, f"{spec}_{family}_sublineages.tsv"), sep="\t", index=False)
    print("Wrote sublineages:", os.path.join(outdir, f"{spec}_{family}_sublineages.tsv"))

    # 9) annotate centroid sizes and compute mean distance zscore for HGT detection
    # collect sizes
    centroid_size = {row.centroid: row.size for row in centroids_df.itertuples()}
    mean_dist = np.mean(distmat, axis=1)
    z = (mean_dist - np.mean(mean_dist)) / (np.std(mean_dist) if np.std(mean_dist)>0 else 1)
    hgt_candidates=[]
    for i,name in enumerate(names):
        size = centroid_size.get(name, 0)
        if z[i] > 2.0 and size >= 3:
            hgt_candidates.append({"centroid":name, "mean_dist":float(mean_dist[i]), "zscore":float(z[i]), "size":size})
    pd.DataFrame(hgt_candidates).to_csv(os.path.join(outdir, f"{spec}_{family}_HGT_candidates.tsv"), sep="\t", index=False)
    print("Wrote HGT candidates file (if any)")

    # 10) FastTree for tree visualization (optional)
    ft = shutil_which("FastTree")
    if ft:
        tree_out = os.path.join(outdir, f"{spec}_{family}_centroids.tree")
        
        # protein alignment: use -wag, no -nt
        cmd = [ft, "-wag", aln_file]
        
        # run with context manager for stdout
        with open(tree_out, "w") as out_f:
            run_cmd(cmd, stdout=out_f)
        
        print("FastTree produced:", tree_out)
    else:
        print("FastTree not found; skip tree generation.")


    print("Pipeline complete. Results in:", outdir)

# small helper
def shutil_which(name):
    from shutil import which
    return which(name)

if __name__ == "__main__":
    main()
