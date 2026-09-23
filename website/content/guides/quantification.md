---
title: Abundance quantification
description: Turn per-read classifications into a taxonomic profile with centrifuger-quant, and export it in Centrifuge, MetaPhlAn, CAMI or Kraken format.
---

`centrifuger` answers "where did each read come from?". `centrifuger-quant` answers the question most analyses actually ask: "what is in this sample, and in what proportion?" It reads a classification file and aggregates it into a taxonomic profile.

```bash
centrifuger-quant -x cfr_idx -c classification.tsv > report.tsv
```

## Inputs

Quantification needs the classification file (`-c`) plus the taxonomy that the classification refers to. The easy way to supply the taxonomy is to point at the same index that produced the classification:

```bash
centrifuger-quant -x cfr_idx -c classification.tsv > report.tsv
```

If the index is not available — say the classification was produced elsewhere, or you would rather not load a 200 GB index just to summarise a TSV — supply the taxonomy files directly instead of `-x`:

```bash
centrifuger-quant \
  --taxonomy-tree taxonomy/nodes.dmp \
  --name-table taxonomy/names.dmp \
  --size-table genome_sizes.tsv \
  -c classification.tsv > report.tsv
```

- `--taxonomy-tree FILE` — NCBI `nodes.dmp`.
- `--name-table FILE` — NCBI `names.dmp`.
- `--size-table FILE` — optional table of contig or genome sizes, used to normalise abundance by genome length.

:::caution Use the matching taxonomy
The taxonomy must be the one the index was built with. A newer `nodes.dmp` may have renamed, merged or deleted nodes that the classification file still refers to.
:::

## The report

The output has seven columns — name, taxID, taxRank, genomeSize, numReads, numUniqueReads and abundance — described in full under [output formats](/reference/output-formats/#quantification-output).

Two of them are easy to confuse:

- **numReads** counts every read assigned anywhere under this taxonomy node. Reads with multiple classifications are distributed evenly among the taxa they were assigned to, so this column can be fractional in spirit even when printed as a count.
- **numUniqueReads** counts only reads assigned unambiguously to a sequence under this node. It is the more conservative evidence, and a large gap between the two columns is a signal that the taxon is hard to distinguish from its neighbours.

**abundance** is the proportion of the sample attributed to this taxon after normalising by genome length, which is what you want when comparing a small virus with a large bacterium.

## Filtering before aggregation

Two options drop weak classifications rather than letting them into the profile:

- `--min-score INT` — ignore reads whose classification score is below the threshold.
- `--min-length INT` — ignore reads whose classified length is below the threshold.

```bash
centrifuger-quant -x cfr_idx -c classification.tsv \
  --min-score 300 --min-length 60 > report.tsv
```

There is no universal cutoff: the right values depend on read length, database breadth and how much you prefer precision to sensitivity. A practical approach is to run without filters first, look at the score distribution in the classification file, and set a threshold that removes the low-score tail.

:::tip Filter at quantification, not classification
Because the filters operate on the classification file, you can try several thresholds without re-running `centrifuger` — which is the expensive step.
:::

## Report formats

`--output-format INT` selects the layout, so the profile can be handed straight to downstream tooling:

| Value | Format |
|-------|--------|
| `0` | Centrifuge report (default) |
| `1` | MetaPhlAn |
| `2` | CAMI |
| `3` | Kraken report |

```bash
# MetaPhlAn-style profile
centrifuger-quant -x cfr_idx -c classification.tsv --output-format 1 > profile.txt

# Kraken-style report, for tools such as Krona or Pavian
centrifuger-quant -x cfr_idx -c classification.tsv --output-format 3 > kreport.txt
```

Since the classification file is the single source of truth, you can emit several formats from one run of `centrifuger` at no extra cost.

## How -k affects quantification

The `-k` option of `centrifuger` — how many distinct assignments are reported per read — carries straight through to the profile. With the default `-k 1`, a read that matches several genomes equally well is reported once; with a larger `-k`, all of its near-equal assignments appear, and `centrifuger-quant` spreads the read across them.

Increasing `-k` therefore gives an **ambiguous but more specific** classification, and it can improve the quantification result, because a read shared between two strains contributes to both rather than being forced onto one. The cost is a larger classification file and more taxa in the report with small fractional support.

If species-level proportions matter to your analysis, it is worth classifying once with `-k 1` and once with, say, `-k 5`, and comparing the resulting profiles.

## A two-step pipeline

```bash
#!/usr/bin/env bash
set -euo pipefail

IDX=/data/indexes/cfr_hpv+gbsarscov2

centrifuger -x "$IDX" -1 sample_1.fq.gz -2 sample_2.fq.gz -t 16 -k 5 \
  > classification.tsv

# default report
centrifuger-quant -x "$IDX" -c classification.tsv > report.tsv

# Kraken-style report for visualisation, same classification file
centrifuger-quant -x "$IDX" -c classification.tsv --output-format 3 > kreport.txt
```

## Next steps

- [Output formats](/reference/output-formats/) — column-by-column reference for both files.
- [FAQ and troubleshooting](/reference/faq/) — interpreting scores, unexpected taxa, empty reports.
