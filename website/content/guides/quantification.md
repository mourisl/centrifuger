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

If the classification used a protein index, also add `--protein`: hit lengths are then in amino acids, and without it every read looks mostly uncovered. With `-x` the index tells `centrifuger-quant` this by itself. Centrifuger reports the hit length in amino acids but the read length in nucleotides; read as nucleotide, every read looks two thirds uncovered and gets the smallest weight, which flattens the abundance estimate. `centrifuger-quant` warns when no hit covers more than a third of its read.

![Abundance of the two example strains with a protein index read as nucleotide or as protein](/assets/protein_weights.png)

The example reads classified against the example index marked as protein, with hit lengths divided by 3. Read as nucleotide (Centrifuger v1.2.0, or taxonomy files without `--protein`), the two strains come out near 50:50. Read as protein (a protein index given with `-x`, or `--protein`), the estimate is close to the one from the nucleotide index.

:::caution Use the matching taxonomy
The taxonomy must be the one the index was built with. A newer `nodes.dmp` may have renamed, merged or deleted nodes that the classification file still refers to.
:::

## The report

The output has nine columns — name, taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance, emNumReads and emNumReadsWeighted — described in full under [output formats](/reference/output-formats/#quantification-output).

Two of them are easy to confuse:

- **numReads** counts every read assigned anywhere under this taxonomy node. Reads with multiple classifications are distributed evenly among the taxa they were assigned to, so this column can be fractional in spirit even when printed as a count.
- **numUniqueReads** counts only reads assigned unambiguously to a sequence under this node. It is the more conservative evidence, and a large gap between the two columns is a signal that the taxon is hard to distinguish from its neighbours.

**abundance** is the proportion of the sample attributed to this taxon after normalising by genome length, which is what you want when comparing a small virus with a large bacterium.

**emNumReads** shares each multi-classified read in proportion to the abundance estimate instead of evenly. When several close relatives share reads, it gives the reads to the taxa the abundance estimate supports, so it is the read count to use alongside **abundance**. Add `--em-reads` to carry it into the MetaPhlAn, CAMI and Kraken report formats.

| Column | What it holds | Use it for |
|--------|---------------|------------|
| numReads | reads under the node, each multi-classified read split evenly; printed as a whole number (the fraction is dropped) | a quick look; comparing with earlier versions |
| numUniqueReads | reads with a single classification under the node | evidence that a taxon is present |
| emNumReads | reads under the node, each multi-classified read shared by the EM abundance estimate; fractional, sums to the classified reads at the root | read counts that agree with the abundance estimate |
| emNumReadsWeighted | the EM's internal read weight | diagnosing the EM only; it is not a read count |
| abundance | proportion of genomes, normalised by genome size | relative abundance (not read counts) |

![One read shared by three strains, as counted by each output](/assets/shared_reads.png)

The figure follows one read that fits three sibling strains (genome sizes 1, 1 and 2 Mb). numReads splits it evenly, 1/3 each, but prints whole numbers, so every strain shows 0. emNumReads shares it by the EM estimate, which divides by genome size: the two 1 Mb strains get half a read each and the 2 Mb strain none. Centrifuger v1.2.0's Kraken report rounded each row on its own, so all three strains got 0 and the taxon reads added up to 0 although the sample has one read. The Kraken report now shares whole reads from the top of the tree down by largest remainder, so one strain gets the read and the taxon reads add up to 1, with or without `--em-reads`.

## Genome sizes

The genome sizes that abundance is divided by come from the index; the genomeSize column and `centrifuger-inspect --size-table` show the same values.

- Each sequence counts its length after dust masking, so a little under its FASTA length unless the index was built with `--no-dust`.
- Sequences of one taxonomy ID are added together as one genome only when their IDs look like an accession with a two-letter prefix and an underscore, and their numbers follow on (`NC_000001`, `NC_000002`).
- When a taxonomy ID has several genomes, the longest is used.
- A higher rank gets the mean of the genomes below it.

A segmented virus whose segment IDs do not follow on, or a genome in contigs named otherwise, is therefore sized by its longest piece, and its abundance comes out too high.

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

# the same, with the EM read counts in the read columns
centrifuger-quant -x cfr_idx -c classification.tsv --output-format 3 --em-reads > kreport_em.txt
```

The Kraken report's read counts are whole numbers that add up: each clade count is its taxon count plus its children's clade counts, and the taxon counts add up to the root's, so Pavian shows the root at 100%.

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
