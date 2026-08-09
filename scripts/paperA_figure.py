#!/usr/bin/env python3
"""Figure 1: naive current-taxid tagging manufactures a false 'before' absence for
taxid-renumbering artifacts, making them indistinguishable from genuinely novel taxa until
the tagging is corrected.

The three percentages below are the article's reported results (genuinely novel taxa
correct-classification rate; artifact panel under naive tagging; artifact panel under
period-accurate tagging, from m1_artifact_verify.py's output). They are hardcoded here
rather than re-read from a results file so this script has no dependency beyond matplotlib.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

groups = [
    "Genuinely novel taxa\n(n=17)",
    "Taxid-renumbering artifacts\nnaive current-taxid tag (n=39)",
    "Taxid-renumbering artifacts\nperiod-accurate tag (n=39)",
]
before = [0.00, 0.00, 99.30]
after = [98.88, 99.30, 99.30]

x = np.arange(len(groups))
width = 0.32

fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=300)
b1 = ax.bar(x - width/2, before, width, label="Before transition (MSL39)", color="#4a6fa5")
b2 = ax.bar(x + width/2, after, width, label="After transition (MSL40)", color="#d97b4f")

for bars in (b1, b2):
    for rect in bars:
        h = rect.get_height()
        ax.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9)

ax.set_ylabel("Mean correct species-level classification (%)")
ax.set_ylim(0, 118)
ax.set_xticks(x)
ax.set_xticklabels(groups, fontsize=9)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=2, frameon=False, fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_title("Naive current-taxid tagging manufactures a false 'before' absence", fontsize=11, pad=45)

fig.tight_layout()
out_path = os.environ.get("PAPERA_FIGURE_OUT", "results/paperA_figure1.png")
os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
fig.savefig(out_path, dpi=300)
print(f"saved: {out_path}")
