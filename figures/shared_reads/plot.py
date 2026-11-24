"""Figure: one read that fits three sibling strains, as counted by each centrifuger-quant output.
Input counts.tsv (from figures/make_data.sh); output shared_reads.png and shared_reads.pdf next to this script."""
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
f_label, f_tick, f_note, f_small = (fm.FontProperties(family=fam, size=s) for s in (11, 10, 9, 7))
colors = ["#77AEB3", "#E5885D", "#C7C4B5", "#A1C2DE", "#B4944B"]

rows = list(csv.DictReader(open(os.path.join(here, "counts.tsv")), delimiter="\t"))
groups = ["Strain_A_1Mb", "Strain_B_1Mb", "Strain_C_2Mb", "taxon_reads_sum"]
group_labels = ["Strain A\n(1 Mb)", "Strain B\n(1 Mb)", "Strain C\n(2 Mb)", "Sum of\ntaxon reads"]
labels = {
    "numReads": "numReads, format 0 (even split 1/3, printed whole)",
    "emNumReads": "emNumReads, format 0 (EM share)",
    "kreport_v1.2.0": "kraken-report, v1.2.0 (each row rounded)",
    "kreport": "kraken-report, this version",
    "kreport_em_reads": "kraken-report --em-reads",
}

fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.set_position([0.09, 0.2, 0.52, 0.7])
n, width = len(rows), 0.16
for i, r in enumerate(rows):
    xs = [g + (i - (n - 1) / 2) * width for g in range(len(groups))]
    ys = [float(r[g]) for g in groups]
    ax.bar(xs, ys, width=width * 0.92, color=colors[i], label=labels[r["output"]], zorder=2)
    for x, y in zip(xs, ys):
        ax.text(x, y + 0.02, f"{y:g}", ha="center", va="bottom", fontproperties=f_small)
ax.axhline(1, color="#555555", linewidth=0.8, linestyle=(0, (3, 3)), zorder=1)
ax.text(len(groups) - 0.5, 1.12, "reads in the sample = 1", ha="right", va="bottom", fontproperties=f_note,
        color="#555555")
ax.set_xticks(range(len(groups)))
ax.set_xticklabels(group_labels, fontproperties=f_tick)
ax.set_ylim(0, 1.3)
ax.set_xlim(-0.6, len(groups) - 0.4)
ax.set_yticks([0, 0.5, 1])
ax.set_yticklabels(["0", "0.5", "1"], fontproperties=f_tick)
ax.set_ylabel("Reads counted", fontproperties=f_label)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, prop=f_note, loc="center left", bbox_to_anchor=(1.02, 0.5))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(here, f"shared_reads.{ext}"), dpi=200, bbox_inches="tight")
