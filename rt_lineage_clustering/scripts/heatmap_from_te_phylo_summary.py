import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# Load dataset
df = pd.read_csv('TE_family_phylo_summary.csv')

# Select metrics to include in the heatmap
metrics = [
    "Lsyl_MPD_z", "Jeff_MPD_z",
    "Lsyl_MNTD_z", "Jeff_MNTD_z",
    "Lsyl_clusters", "Jeff_clusters",
    "between"
]

# Subset and set index
heatmap_df = df.set_index("family")[metrics]

# Scale values (z-scoring)
scaled = pd.DataFrame(
    StandardScaler().fit_transform(heatmap_df),
    index=heatmap_df.index,
    columns=heatmap_df.columns
)

# Create a clustered heatmap
plt.figure(figsize=(12, 8))
sns.clustermap(
    scaled,
    cmap="coolwarm",
    linewidths=0.5,
    figsize=(14,10),
    col_cluster=True,
    row_cluster=True,
    dendrogram_ratio=(0.2, 0.2),
    cbar_pos=(0.04, 0.8, 0.02, 0.18)
)

plt.suptitle("TE Family Dynamics Heatmap (scaled values)", y=1.02, fontsize=16)

# Save to file
output_path = "TE_dynamics_heatmap.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')

output_path
