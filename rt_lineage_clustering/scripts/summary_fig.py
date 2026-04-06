import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import sys

csv = sys.argv[1]

# Load file
df = pd.read_csv(csv)

plt.figure(figsize=(10,8))

# Scatter
plt.scatter(df['Lsyl_MPD_z'], df['Jeff_MPD_z'], s=60)

# Collect labels
texts = []
for i, row in df.iterrows():
    texts.append(
        plt.text(
            row['Lsyl_MPD_z'],
            row['Jeff_MPD_z'],
            row['family'],
            fontsize=9
        )
    )

# Adjust labels to prevent overlap
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', lw=0.5, color='gray'),
    force_points=0.4,
    force_text=0.4,
    expand_points=(1.2, 1.4),
    expand_text=(1.2, 1.4),
)

plt.xlabel('Lsyl MPD z-score')
plt.ylabel('Jeff MPD z-score')
plt.title('MPD z-score comparison between Lsyl and Jeff')
plt.axvline(0,linewidth=1,color ="orange",)
plt.axhline(0,linewidth=1,color="orange")
plt.grid(True, alpha=0.3)

plt.savefig('TE_MPD_z_scores.png', dpi=300, bbox_inches='tight')

