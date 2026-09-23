---
title: Output formats
description: Column-by-column description of the Centrifuger classification TSV and the centrifuger-quant abundance report.
---

Centrifuger produces two tab-separated files: a per-read classification from `centrifuger`, and a per-taxon profile from `centrifuger-quant`. Both go to standard output.

## Inputs

The primary input to Centrifuger is the index of the genome or protein database (`-x`) together with gzipped or uncompressed FASTQ read files — `-1`/`-2` for paired-end data, `-u` for single-end, `-i` for interleaved.

## Classification output

`centrifuger` writes one row per read assignment, with eight columns:

```text
readID    seqID          taxID     score   2ndBestScore   hitLength   queryLength   numMatches
1_1       MT019531.1     2697049   4225    0              80          80            1
```

| # | Column | Description |
|---|--------|-------------|
| 1 | `readID` | The read ID from the raw sequencing read, e.g. `1_1`. |
| 2 | `seqID` | The sequence ID of the genomic sequence the read was classified to, e.g. `MT019531.1`. |
| 3 | `taxID` | The taxonomic ID of the sequence in column 2, e.g. `2697049`. |
| 4 | `score` | The score for the classification: the weighted sum of hits, e.g. `4225`. |
| 5 | `2ndBestScore` | The score of the next best classification, e.g. `0`. |
| 6 | `hitLength` | The number of base pairs of the read that match the genomic sequence found by Centrifuger, e.g. `80`. |
| 7 | `queryLength` | The length of the read, or the combined length of the mate pair, e.g. `80`. |
| 8 | `numMatches` | The number of classifications made for this read, i.e. how many assignments appear in the output, e.g. `1`. |

### Reading the file

- **A read can appear more than once.** With `-k` greater than 1, near-equal assignments each get their own row; `numMatches` tells you how many rows to expect for that read.
- **`score` versus `2ndBestScore`.** A high score with a `2ndBestScore` of 0 is an unambiguous assignment. Two nearly equal scores mean the read cannot distinguish the two references — the usual situation for conserved regions and closely related strains.
- **`hitLength` versus `queryLength`.** Their ratio is how much of the read actually matched. A long read with a short hit length is weak evidence even when the score looks large.
- **Unclassified reads** are not represented by a classification row. To collect them, pass `--un` to `centrifuger`.

When barcode and UMI extraction are enabled, those values ride along with the read so that assignments can be grouped by cell afterwards — see [single-cell and barcoded data](/guides/single-cell/).

### Quick summaries

```bash
# how many reads were classified
cut -f1 classification.tsv | tail -n +2 | sort -u | wc -l

# most frequently hit taxa
tail -n +2 classification.tsv | cut -f3 | sort | uniq -c | sort -rn | head

# keep only confident, unambiguous assignments
awk -F'\t' 'NR==1 || ($4 > 300 && $5 == 0)' classification.tsv > confident.tsv
```

## Quantification output

`centrifuger-quant` estimates the abundance of each taxonomy ID. The report has seven columns:

```text
name                                                          taxID    taxRank  genomeSize  numReads  numUniqueReads  abundance
Legionella_pneumophila_subsp._pneumophila_str._Philadelphia_1  272624   strain   3397753     50        48              0.392641
```

| # | Column | Description |
|---|--------|-------------|
| 1 | `name` | The name of a genome, or the name corresponding to the taxonomic ID in column 2 at a rank higher than strain. |
| 2 | `taxID` | The taxonomic ID. |
| 3 | `taxRank` | The taxonomic rank, e.g. `strain`. |
| 4 | `genomeSize` | The length of the genome sequence. |
| 5 | `numReads` | The number of reads classified to some genomic sequence under this taxonomy node. Multi-classified reads are evenly distributed. |
| 6 | `numUniqueReads` | The number of reads uniquely classified to a genomic sequence under this taxonomy node. |
| 7 | `abundance` | The proportion of this genome, normalised by its genomic length. |

### Reading the report

- **`numReads` against `numUniqueReads`.** A taxon with many reads but few unique ones is hard to distinguish from its neighbours; one with a high unique count is well supported.
- **`abundance` is length-normalised**, so it is comparable between a small virus and a large bacterium in a way that raw read counts are not.
- **Rows exist at several ranks.** Reads that cannot be resolved to strain are attributed to the lowest rank that does explain them, so a profile mixes strain, species and higher-rank rows.

### Alternative formats

`--output-format` re-lays the same information for downstream tools:

| Value | Format | Typical consumer |
|-------|--------|------------------|
| `0` | Centrifuge report (default) | the columns described above |
| `1` | MetaPhlAn | MetaPhlAn-compatible profile tooling |
| `2` | CAMI | CAMI benchmarking and profile comparison |
| `3` | Kraken report | Krona, Pavian, Bracken-style tooling |

```bash
centrifuger-quant -x cfr_idx -c classification.tsv --output-format 3 > kreport.txt
```

Since these are all derived from the classification file, emitting several costs nothing beyond a second `centrifuger-quant` run.

## Index files

`centrifuger-build` writes `<prefix>.1.cfr` through `<prefix>.4.cfr`. These four files are the whole index; the reference FASTA and intermediate build files can be deleted once the build succeeds. Pass the prefix — the path without `.1.cfr` — to `-x`.

With `--checkpoint`, the build also writes `<prefix>_checkpoint.[123]`, which exist only to resume an interrupted build and can be removed afterwards.
