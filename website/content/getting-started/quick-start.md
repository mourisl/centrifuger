---
title: Quick start
description: Classify your first sample end to end, using either the bundled example data or a pre-built index.
---

This page takes you from a fresh installation to an abundance profile. There are two routes: a tiny bundled example that runs in seconds and proves the installation works, and a real run against a pre-built database.

:::note Before you begin
Centrifuger must be installed and on your `PATH`. See [installation](/getting-started/installation/) if it is not.
:::

![Centrifuger workflow: build an index once, classify each sample, then quantify](/assets/workflow.png)

The index carries the taxonomy, names and genome sizes, so classification and quantification need only `-x`. Quantification takes seconds, so it can be rerun with other options without classifying again.

## Run the bundled example

The source distribution contains an `example/` directory with a small reference FASTA, an NCBI-style taxonomy, and a pair of FASTQ files. Move into that directory first:

```bash
cd /path/to/centrifuger/example
```

### 1. Build a small index

```bash
../centrifuger-build -r ref.fa \
  --taxonomy-tree nodes.dmp \
  --name-table names.dmp \
  --conversion-table ref_seqid.map \
  -o cfr_ref_idx
```

Afterwards you should see the index files `cfr_ref_idx.*.cfr` in the folder.

### 2. Classify the reads

```bash
../centrifuger -1 example_1.fq -2 example_2.fq -x cfr_ref_idx > output.tsv
```

`output.tsv` should closely match the `example_class.out` file shipped alongside it. If it does, your installation is working. All 100 read pairs are classified.

### 3. Summarise into an abundance profile

```bash
../centrifuger-quant -x cfr_ref_idx -c output.tsv > report.tsv
```

The example holds two strains of *Legionella pneumophila*. Their rows in `report.tsv` (columns name and taxRank left out):

```
taxID   genomeSize  numReads  numUniqueReads  abundance  emNumReads  emNumReadsWeighted
272624  3353119     47        47              0.3478411  47.000      6.607
297246  3458445     47        47              0.6521589  47.000      12.722
```

With the default `-k 1`, the six reads that fit both strains equally well are assigned to the species above them, so each strain has 47 reads, all unique. Classify with `-k 5` instead and those six reads keep both assignments: numReads then splits them evenly (50 and 50) and emNumReads shares them by the abundance estimate (49.089 and 50.911); see [abundance quantification](/guides/quantification/).

For a Kraken-style report whose read columns are the EM counts:

```bash
../centrifuger-quant -x cfr_ref_idx -c output.tsv --output-format 3 --em-reads > kreport.txt
```

## Classify a real sample

For real data you need a database that covers the organisms you care about. Downloading a published index is far quicker than building one.

### 1. Get an index

```bash
centrifuger-download cfr_hpv+gbsarscov2
```

`cfr_hpv+gbsarscov2` is RefSeq human, bacteria, archaea and virus plus SARS-CoV-2 variants from GenBank — about 41 GB. Other databases, including GTDB r232, NCBI `core_nt` and NCBI `nr`, are listed on the [pre-built indexes](/guides/prebuilt-indexes/) page.

:::caution Check your memory first
The index is loaded into RAM. Make sure the machine has at least as much memory as the index size listed in the download table, plus a little headroom.
:::

### 2. Classify

Paired-end reads, using eight threads:

```bash
centrifuger -x cfr_hpv+gbsarscov2 \
  -1 sample_1.fq.gz -2 sample_2.fq.gz \
  -t 8 > classification.tsv
```

Single-end data uses `-u` instead. Long reads are normally single-end, so they take this form too:

```bash
centrifuger -x cfr_hpv+gbsarscov2 -u nanopore.fq.gz -t 8 > classification.tsv
```

Each row of `classification.tsv` is one read assignment. The columns are described in [output formats](/reference/output-formats/#classification-output).

### 3. Quantify

```bash
centrifuger-quant -x cfr_hpv+gbsarscov2 -c classification.tsv > report.tsv
```

`report.tsv` holds one row per taxon with read counts and a length-normalised abundance. To hand the profile to a downstream tool, ask for a different format — for example a Kraken-style report:

```bash
centrifuger-quant -x cfr_hpv+gbsarscov2 -c classification.tsv \
  --output-format 3 > kraken_style_report.tsv
```

## A complete, copy-pasteable pipeline

```bash
#!/usr/bin/env bash
set -euo pipefail

IDX=cfr_hpv+gbsarscov2
THREADS=8

centrifuger -x "$IDX" -1 sample_1.fq.gz -2 sample_2.fq.gz \
  -t "$THREADS" --un unclassified \
  > classification.tsv

centrifuger-quant -x "$IDX" -c classification.tsv \
  > report.tsv
```

`--un unclassified` additionally writes the reads that could not be assigned to `unclassified_1.fq.gz` and `unclassified_2.fq.gz`, which is useful when you suspect the database is missing something.

## Where to go next

- [Building an index](/guides/building-an-index/) — when no published database fits your project.
- [Classifying reads](/guides/classification/) — sample sheets, wildcards, `-k`, low-complexity masking and the rest.
- [Abundance quantification](/guides/quantification/) — filtering and report formats.
- [Single-cell and barcoded data](/guides/single-cell/) — 10x Genomics, UMIs and barcode translation.
