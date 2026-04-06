import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import sys

# =====================
# Styling (global)
# =====================
plt.rcParams.update({
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})

# =====================
# Input
# =====================
csv = sys.argv[1]

# Load file
df = pd.read_csv(csv)

# Filter only families with significant permutation differences
sig_thresh = 0.05
df_sig = df[df['perm_p'] < sig_thresh]

if df_sig.empty:
    print(f"No families pass the significance threshold (perm_p < {sig_thresh}).")
    sys.exit()

# =====================
# MPD plot
# =====================
plt.figure(figsize=(5, 5))

colors = []
for x, y in zip(df_sig['Lsyl_MPD_z'], df_sig['Jeff_MPD_z']):
    if x < 0 and y < 0:
        colors.append('red')    # burst in both
    elif x < 0 and y > 0:
        colors.append('blue')   # Lsyl burst only
    elif x > 0 and y < 0:
        colors.append('green')  # Jeff burst only
    else:
        colors.append('grey')   # old/suppressed

plt.scatter(
    df_sig['Lsyl_MPD_z'],
    df_sig['Jeff_MPD_z'],
    c=colors,
    s=70
)

# Add text labels
texts = []
for _, row in df_sig.iterrows():
    texts.append(
        plt.text(
            row['Lsyl_MPD_z'],
            row['Jeff_MPD_z'],
            row['family'],
            fontsize=12,
            fontweight='bold'
        )
    )

# Adjust label positions to avoid overlap
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', lw=0.6, color='gray'),
    force_points=0.3,
    force_text=0.4,
    expand_points=(1.2, 1.4),
    expand_text=(1.2, 1.4)
)

# Axes and formatting
plt.xlabel('Lsyl MPD z-score')
plt.ylabel('Jeff MPD z-score')
plt.title(f'MPD z-score comparison — Significant families only (perm_p < {sig_thresh})')

plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.axvline(0, color='black', linestyle='--', linewidth=1)

plt.grid(True, alpha=0.3)

# Save plot
plt.savefig(
    'TE_MPD_z_scores_significant_formatted.png',
    dpi=300,
    bbox_inches='tight'
)
print("Saved: TE_MPD_z_scores_significant_formatted.png")

plt.close()

# =====================
# MNTD plot
# =====================
plt.figure(figsize=(5, 5))

colors = []
for x, y in zip(df_sig['Lsyl_MNTD_z'], df_sig['Jeff_MNTD_z']):
    if x < 0 and y < 0:
        colors.append('red')    # burst in both
    elif x < 0 and y > 0:
        colors.append('blue')   # Lsyl burst only
    elif x > 0 and y < 0:
        colors.append('green')  # Jeff burst only
    else:
        colors.append('grey')   # old/suppressed

plt.scatter(
    df_sig['Lsyl_MNTD_z'],
    df_sig['Jeff_MNTD_z'],
    c=colors,
    s=70
)

# Add text labels
texts = []
for _, row in df_sig.iterrows():
    texts.append(
        plt.text(
            row['Lsyl_MNTD_z'],
            row['Jeff_MNTD_z'],
            row['family'],
            fontsize=12,
            fontweight='bold'
        )
    )

# Adjust label positions to avoid overlap
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', lw=0.6, color='gray'),
    force_points=0.3,
    force_text=0.4,
    expand_points=(1.2, 1.4),
    expand_text=(1.2, 1.4)
)

# Axes and formatting
plt.xlabel('Lsyl MNTD z-score')
plt.ylabel('Jeff MNTD z-score')
plt.title(f'MNTD z-score comparison — Significant families only (perm_p < {sig_thresh})')

plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.axvline(0, color='black', linestyle='--', linewidth=1)

plt.grid(True, alpha=0.3)

# Save plot
plt.savefig(
    'TE_MNTD_z_scores_significant_formatted.png',
    dpi=300,
    bbox_inches='tight'
)
print("Saved: TE_MNTD_z_scores_significant_formatted.png")

plt.close()
