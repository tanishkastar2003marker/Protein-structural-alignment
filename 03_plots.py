
import pandas as pd
import matplotlib.pyplot as plt
import os

pairs = [
    "Lcasei_DHFR_vs_Ecoli_DHFR",
    "Haemoglobin-a_vs_Myoglobin",
    "T4_Lysozyme_vs_Hen_Lysozyme",
    "Cytochrome-c_vs_Cytochrome-c2",
    "Azurin_vs_Plastocyanin"
]

nice_names = [
    "L.casei vs E.coli DHFR",
    "Haemoglobin vs Myoglobin",
    "T4 vs Hen Lysozyme",
    "Cytochrome-c vs c2",
    "Azurin vs Plastocyanin"
]

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()

for idx, (pair, name) in enumerate(zip(pairs, nice_names)):
    filepath = f"results/{pair}.csv"
    if not os.path.exists(filepath):
        print(f"Missing: {filepath}")
        continue

    df = pd.read_csv(filepath)
    ax = axes[idx]

    ax.plot(df['τ (Å)'], df['τ-equivalences'],
            'b-o', linewidth=2, label='τ-equivalences')
    ax.plot(df['τ (Å)'], df['δ-equivalences'],
            'r-s', linewidth=2, label='δ-equivalences')

    ax.set_title(name, fontsize=10, fontweight='bold')
    ax.set_xlabel('τ value (Å)')
    ax.set_ylabel('Number of equivalences')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

# hide the empty 6th subplot
axes[5].set_visible(False)

plt.suptitle(
    'τ and δ Equivalences vs Tolerance Value\n'
    'Reproducing Subbarao & Haneef (1991)',
    fontsize=13, fontweight='bold'
)
plt.tight_layout()
plt.savefig('results/equivalences_plot.png',
            dpi=150, bbox_inches='tight')
print("Plot saved to results/equivalences_plot.png")
plt.show()