"""Figure: the Centrifuger workflow from reference genomes to a taxonomic profile.
No input data; the boxes and arrows are defined below. Output workflow.png and workflow.pdf next to this script."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

matplotlib.rcParams["mathtext.default"] = "regular"
matplotlib.rcParams["pdf.fonttype"] = 42
here = os.path.dirname(os.path.abspath(__file__))
fam = ["Liberation Sans", "DejaVu Sans"]
f_box = fm.FontProperties(family=fam, size=9)
f_prog = fm.FontProperties(family=fam, size=10)
f_note = fm.FontProperties(family=fam, size=8)
f_head = fm.FontProperties(family=fam, size=9)

FILE, PROG, OUT = "#E8E6DD", "#77AEB3", "#F2C9AE"
EDGE = "#555555"

fig = plt.figure(figsize=(9.6, 5.0))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 96)
ax.set_ylim(0, 50)
ax.axis("off")

boxes = {}


def box(key, x, y, w, h, text, color, font=f_box):
    """Rounded box with its lower-left corner at (x, y)."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                facecolor=color, edgecolor=EDGE, linewidth=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontproperties=font, linespacing=1.3)
    boxes[key] = (x, y, w, h)


def side(key, where):
    x, y, w, h = boxes[key]
    return {"l": (x - 0.3, y + h / 2), "r": (x + w + 0.3, y + h / 2),
            "t": (x + w / 2, y + h + 0.3), "b": (x + w / 2, y - 0.3)}[where]


def arrow(a, b, rad=0.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=10, color=EDGE, linewidth=0.9,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


def note(x, y, text):
    ax.text(x, y, text, fontproperties=f_note, va="top", color="#333333", linespacing=1.3)


ax.text(1, 48, "1. Build the index (once per database)", fontproperties=f_head, va="center")
ax.text(1, 21.5, "2. Classify each sample", fontproperties=f_head, va="center")
ax.text(58, 21.5, "3. Quantify (seconds; rerun freely)", fontproperties=f_head,
        va="center")

# Build row: three inputs feed centrifuger-build
box("ref", 1, 39.5, 21, 5.5, "Reference sequences, FASTA\n(-r FILE, or -l FILE list)", FILE)
box("map", 1, 32.5, 21, 5.5, "seqID → taxID table\n(--conversion-table)", FILE)
box("tax", 1, 25.5, 21, 5.5, "nodes.dmp, names.dmp\n(--taxonomy-tree, --name-table)", FILE)
box("build", 28, 32.5, 16, 5.5, "centrifuger-build", PROG, f_prog)
box("idx", 50, 32.5, 15, 5.5, "Index\nprefix.1.cfr … prefix.4.cfr", OUT)
note(67.5, 38.6, "The index carries the sequences,\ntaxonomy, names and genome\nsizes, so later steps need only\n-x prefix (no .dmp files).")
arrow(side("ref", "r"), (28 - 0.3, 36.6), rad=-0.12)
arrow(side("map", "r"), side("build", "l"))
arrow(side("tax", "r"), (28 - 0.3, 33.6), rad=0.12)
arrow(side("build", "r"), side("idx", "l"))

# Classify and quantify row
box("reads", 1, 10, 16, 6, "Reads, FASTQ (.gz)\n-1/-2, -u or -i", FILE)
box("cfr", 21, 10, 14, 6, "centrifuger\n-x prefix -k 5", PROG, f_prog)
box("cls", 39, 10, 15, 6, "classification.tsv\none row per hit", OUT)
box("quant", 58, 10, 15, 6, "centrifuger-quant\n-x prefix -c", PROG, f_prog)
arrow(side("reads", "r"), side("cfr", "l"))
arrow(side("cfr", "r"), side("cls", "l"))
arrow(side("cls", "r"), side("quant", "l"))
arrow(side("idx", "b"), side("cfr", "t"), rad=-0.12)
arrow(side("idx", "b"), side("quant", "t"))

# Outputs
box("f0", 77, 17.5, 18.5, 6.5, "report.tsv (default)\nnumReads, numUniqueReads,\nabundance, emNumReads", OUT)
box("f12", 77, 9.75, 18.5, 5.5, "--output-format 1 or 2\nMetaPhlAn or CAMI profile", OUT)
box("f3", 77, 1, 18.5, 6.5, "--output-format 3\nKraken report (Pavian, Krona)\nadd --em-reads for EM counts", OUT)
arrow(side("quant", "r"), side("f0", "l"), rad=0.15)
arrow(side("quant", "r"), side("f12", "l"))
arrow(side("quant", "r"), side("f3", "l"), rad=-0.15)

note(21, 7.2, "-k 5: report up to five best assignments per read;\ncentrifuger-quant shares such reads between taxa.")
note(21, 3.2, "Protein index: build with --protein; centrifuger-quant\nreads the index type from -x (else pass --protein).")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(here, f"workflow.{ext}"), dpi=200)
