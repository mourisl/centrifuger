"""Figure: abundance of the two example strains when a protein index's hits are read as nucleotide or as protein.
Input weights.tsv (from figures/make_data.sh); output protein_weights.png and protein_weights.pdf next to this script."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

matplotlib.rcParams["mathtext.default"] = "regular"
matplotlib.rcParams["pdf.fonttype"] = 42
here = os.path.dirname(os.path.abspath(__file__))
fam = ["Liberation Sans", "DejaVu Sans"]
f_label, f_tick, f_note = (fm.FontProperties(family=fam, size=s) for s in (11, 10, 9))
colors = ["#77AEB3", "#E5885D"]

rows = list(csv.DictReader(open(os.path.join(here, "weights.tsv")), delimiter="\t"))
runs = ["nt", "aa_as_nt", "aa_protein"]
run_labels = ["Nucleotide index\n(reference)", "Protein index read\nas nucleotide\n(v1.2.0; no --protein)",
              "Protein index\n(-x protein index\nor --protein)"]
strains = [("272624", "Strain Philadelphia 1"), ("297246", "Strain Paris")]
abund = {(r["run"], r["taxID"]): float(r["abundance"]) for r in rows}

fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.set_position([0.11, 0.24, 0.58, 0.66])
width = 0.34
for i, (tid, name) in enumerate(strains):
    xs = [k + (i - 0.5) * width for k in range(len(runs))]
    ys = [abund[(run, tid)] for run in runs]
    ax.bar(xs, ys, width=width * 0.92, color=colors[i], label=name, zorder=2)
    for x, y in zip(xs, ys):
        ax.text(x, y + 0.015, f"{y:.2f}", ha="center", va="bottom", fontproperties=f_note)
ax.set_xticks(range(len(runs)))
ax.set_xticklabels(run_labels, fontproperties=f_tick)
ax.set_xlim(-0.6, len(runs) - 0.4)
ax.set_ylim(0, 0.8)
ax.set_yticks([0, 0.25, 0.5, 0.75])
ax.set_yticklabels(["0", "0.25", "0.5", "0.75"], fontproperties=f_tick)
ax.set_ylabel("Abundance", fontproperties=f_label)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, prop=f_note, loc="center left", bbox_to_anchor=(1.02, 0.5))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(here, f"protein_weights.{ext}"), dpi=200, bbox_inches="tight")
