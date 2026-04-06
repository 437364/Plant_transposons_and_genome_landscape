import re
import glob
import pandas as pd

def extract_float(pattern, text):
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None

def extract_two(pattern, text):
    m = re.search(pattern, text)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None

def parse_within(text):
    # capture two numbers after species names
    m = re.search(r"Within-species mean patristic distances:\s+Lsyl\s+Jeff\s+([\d\.eE\-]+)\s+([\d\.eE\-]+)", text)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None

def parse_mpd_block(text, species):
    # species lines look like:
    # Lsyl   330 0.1532368 0.2511729 ... -18.898868
    pat = rf"{species}\s+\d+\s+([\d\.eE\-]+)\s+[\d\.eE\-]+\s+[\d\.eE\-]+\s+\d+\s+([\d\.eE\-]+)"
    m = re.search(pat, text)
    if not m:
        return None, None, None
    obs = float(m.group(1))
    z = float(m.group(2))

    # get p from "mpd.obs.p" table
    pat2 = rf"{species}\s+([\d\.eE\-]+)\s+\d+"
    m2 = re.search(pat2, text[text.find("mpd.obs.p"):])
    p = float(m2.group(1)) if m2 else None

    return obs, z, p

def parse_mntd_block(text, species):
    pat = rf"{species}\s+\d+\s+([\d\.eE\-]+)\s+[\d\.eE\-]+\s+[\d\.eE\-]+\s+\d+\s+([\d\.eE\-]+)"
    m = re.search(pat, text[text.find("MNTD / NTI"):])
    if not m:
        return None, None, None
    obs = float(m.group(1))
    z = float(m.group(2))

    # p-value
    pat2 = rf"{species}\s+([\d\.eE\-]+)\s+\d+"
    m2 = re.search(pat2, text[text.find("mntd.obs.p"):])
    p = float(m2.group(1)) if m2 else None

    return obs, z, p

def parse_clusters(text):
    # find cluster block
    if "Clusters (founder groups)" not in text:
        return None, None

    clust_text = text[text.find("Clusters (founder groups)"):]
    
    # Each cluster line: "   1     0  229"
    lines = re.findall(r"^\s*\d+\s+\d+\s+\d+", clust_text, flags=re.M)

    Jeff_c = 0
    Lsyl_c = 0

    for ln in lines:
        parts = ln.split()
        _, jeff, lsyl = parts
        if int(jeff) > 0:
            Jeff_c += 1
        if int(lsyl) > 0:
            Lsyl_c += 1

    return Jeff_c, Lsyl_c


# ==========================
# MAIN
# ==========================

records = []

for file in glob.glob("*_phylo_analysis.txt"):
    with open(file) as f:
        text = f.read()

    family = file.replace("_TE_phylo_analysis.txt", "")

    n_tips = extract_float(r"Loaded tree with (\d+) tips", text)
    between = extract_float(r"Between-species mean patristic distance:\s*([\d\.eE\-]+)", text)

    Lsyl_within, Jeff_within = parse_within(text)

    # MPD
    Lsyl_mpd_obs, Lsyl_mpd_z, Lsyl_mpd_p = parse_mpd_block(text, "Lsyl")
    Jeff_mpd_obs, Jeff_mpd_z, Jeff_mpd_p = parse_mpd_block(text, "Jeff")

    # MNTD
    Lsyl_mntd_obs, Lsyl_mntd_z, Lsyl_mntd_p = parse_mntd_block(text, "Lsyl")
    Jeff_mntd_obs, Jeff_mntd_z, Jeff_mntd_p = parse_mntd_block(text, "Jeff")

    perm_p = extract_float(r"Permutation P-value:\s*([\d\.eE\-]+)", text)

    # clusters
    Jeff_clusters, Lsyl_clusters = parse_clusters(text)

    records.append(dict(
        family=family,
        n_tips=n_tips,
        Lsyl_within=Lsyl_within,
        Jeff_within=Jeff_within,
        between=between,

        Lsyl_MPD_z=Lsyl_mpd_z,
        Jeff_MPD_z=Jeff_mpd_z,
        Lsyl_MPD_p=Lsyl_mpd_p,
        Jeff_MPD_p=Jeff_mpd_p,

        Lsyl_MNTD_z=Lsyl_mntd_z,
        Jeff_MNTD_z=Jeff_mntd_z,
        Lsyl_MNTD_p=Lsyl_mntd_p,
        Jeff_MNTD_p=Jeff_mntd_p,

        perm_p=perm_p,

        Jeff_clusters=Jeff_clusters,
        Lsyl_clusters=Lsyl_clusters
    ))

df = pd.DataFrame(records)
df.to_csv("TE_family_phylo_summary.csv", index=False)

print("\nWritten TE_family_phylo_summary.csv")
print(df)
