---
title: Command-line interface
description: Complete option reference for centrifuger-build, centrifuger, centrifuger-quant and centrifuger-download.
---

Every Centrifuger program prints its usage message when run without arguments. This page collects those options in one place.

## centrifuger-build

Builds a searchable index from reference sequences and an NCBI-style taxonomy. See [building an index](/guides/building-an-index/) for worked examples.

```text
Usage: ./centrifuger-build [OPTIONS]
```

### Required

| Option | Description |
|--------|-------------|
| `-r FILE` | Reference sequence file. Repeat `-r` to give more than one input file. |
| `-l FILE` | List of reference sequence files, one path per row. A second column may carry the taxonomy ID mapping. Use instead of `-r`. |
| `--taxonomy-tree FILE` | Taxonomy tree, i.e. the `nodes.dmp` file. |
| `--name-table FILE` | Name table, i.e. the `names.dmp` file. |

`-r` and `-l` are alternatives; supply one of them.

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `--conversion-table FILE` | — | seqID to taxID conversion file. When not set, `-l` is expected and its file must have two columns, `file taxID`. |
| `-o STRING` | `centrifuger` | Output prefix. |
| `-t INT` | `1` | Number of threads. |
| `--protein` | genome | Reference consists of protein sequences. |
| `--build-mem STR` | not used | Automatically infer `--bmax` and `--dcv` to match a memory constraint. Accepts `T`, `G`, `M`, `K` suffixes. |
| `--bmax INT` | `16777216` | Block size for blockwise suffix array sorting. |
| `--dcv INT` | `4096` | Difference cover period. |
| `--offrate INT` | `4` | SA/offset is sampled every 2^INT BWT characters. |
| `--subset-tax INT` | `0` | Only consider the subset of input genomes under this taxonomy node. |
| `--concat-tax-genome` | not used | Concatenate genomes sharing a taxID and discard the seqID information. |
| `--ignore-uncategorized-genome` | include all | Ignore genomes whose seqID or taxID is missing or uncategorised. |
| `--checkpoint` | not used | Write checkpoint files (`[output_prefix]_checkpoint.[123]`) so an interrupted build can resume. |

:::caution Defaults are tuned for small databases
The default `--bmax` and `--dcv` may be inefficient for larger genome databases. Use `--build-mem` with a rough estimate of the available memory instead of tuning them by hand.
:::

## centrifuger

Classifies reads against an index. Output is written to standard output. See [classifying reads](/guides/classification/).

```text
Usage: ./centrifuger [OPTIONS] > classification.tsv
```

### Required

| Option | Description |
|--------|-------------|
| `-x FILE` | Index prefix. |
| `-1 FILE -2 FILE` | Paired-end read files. |
| `-u FILE` | Single-end read file. |
| `-i FILE` | Interleaved paired-end read file. |
| `--sample-sheet FILE` | List of sample files. Each row: `read1 read2 barcode UMI output`. Use a dot (`.`) where there is no such file. |

`-x` is always required, together with exactly one of the read-input forms.

### Optional

| Option | Default | Description |
|--------|---------|-------------|
| `-t INT` | `1` | Number of threads. |
| `-k INT` | `1` | Report up to this many distinct, primary assignments for each read pair. |
| `--un STR` | — | Write unclassified reads to files with this prefix, e.g. `<str>_1/2.fq.gz`. |
| `--cl STR` | — | Write classified reads to files with this prefix. |
| `--barcode STR` | — | Path to the barcode file. |
| `--UMI STR` | — | Path to the UMI file. |
| `--read-format STR` | — | Format for read, barcode and UMI files, e.g. `r1:0:-1,r2:0:-1,bc:0:15,um:16:-1`. |
| `--barcode-whitelist STR` | — | Path to the barcode whitelist file. Requires `--barcode`. |
| `--barcode-translate STR` | — | Path to the barcode translation file. |
| `--min-hitlen INT` | auto | Minimum length of partial hits. |
| `--hitk-factor INT` | `40` | Resolve at most `INT × k` entries for each hit. Use `0` for no restriction. |
| `--consider-secondary INT,FLOAT` | `2000,0.995` | Consider a secondary hit when its hit length ≥ `INT` and its score ≥ `FLOAT × best_score`. |
| `--no-dust` | mask | Do not DUST-mask low-complexity regions of reads. |
| `--merge-readpair` | — | Merge overlapping paired-end reads and trim adapters. |

### The --read-format specification

Each comma-separated field takes the form:

```text
[r1|r2|bc|um]:start:end:strand
```

`start` and `end` are 0-based and inclusive; `-1` means the end of the read. `strand` is `+` or `-`, may be omitted when `+`, and is ignored for `r1` and `r2`. Repeating a field concatenates non-consecutive segments, e.g. `bc:0:15,bc:32:-1`.

To read a barcode or UMI out of the FASTQ header comment instead of the sequence:

```text
[bc|um]:hd:field:start:end:strand
```

`hd` is a keyword selecting the header comment. `field` is either a 0-based field index within the comment (the read ID excluded) or a string prefix to search for. Full explanation and examples in [single-cell and barcoded data](/guides/single-cell/).

## centrifuger-quant

Aggregates a classification file into a taxonomic profile. See [abundance quantification](/guides/quantification/).

```text
Usage: ./centrifuger-quant [OPTIONS] > report.tsv
```

### Required

| Option | Description |
|--------|-------------|
| `-x FILE` | Index prefix. |
| `-c FILE` | Classification file. |

When `-x` is not given, the taxonomy must be supplied directly:

| Option | Description |
|--------|-------------|
| `--taxonomy-tree FILE` | Taxonomy tree, i.e. the `nodes.dmp` file. |
| `--name-table FILE` | Name table, i.e. the `names.dmp` file. |
| `--size-table FILE` | Optional table of contig or genome sizes. |

### Optional

| Option | Description |
|--------|-------------|
| `--min-score INT` | Only consider reads with a score of at least this value. |
| `--min-length INT` | Only consider reads with a classified length of at least this value. |
| `--output-format INT` | Output format: `0` Centrifuge (default), `1` MetaPhlAn, `2` CAMI, `3` Kraken report. |

## centrifuger-download

Downloads reference sequences and taxonomy from NCBI, and fetches published indexes. It follows the conventions of `centrifuge-download`, so [Centrifuge's documentation](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building) applies with the program renamed.

```bash
# fetch a pre-built index by title
centrifuger-download cfr_hpv+gbsarscov2

# NCBI taxonomy into taxonomy/
centrifuger-download -o taxonomy taxonomy

# RefSeq archaea, bacteria and viruses into library/, capturing the seqID→taxID map
centrifuger-download -o library -d "archaea,bacteria,viral" refseq > seqid2taxid.map
```

The options used across this documentation:

| Option | Description |
|--------|-------------|
| `-o DIR` | Output directory. |
| `-d STR` | Comma-separated list of NCBI domains, e.g. `archaea,bacteria,viral`, `vertebrate_mammalian`. |
| `-t INT` | Restrict to this taxonomy ID, e.g. `9606` for human, `10090` for mouse. |
| `-a STR` | Restrict to this assembly level, e.g. `Chromosome`. |
| `-c STR` | Restrict to this RefSeq category, e.g. `reference genome`. |

:::note Full option list
Run `centrifuger-download` with no arguments to print the complete, version-accurate usage message for your installation.
:::

## Exit status and streams

- Classification and quantification results go to **standard output**; progress and diagnostics go to standard error. Redirecting stdout to a file therefore keeps the messages visible.
- A filename of `-` reads from standard input, which is how `samtools` output is piped into `centrifuger`.
- Wildcards in `-1`, `-2` and `-u` are expanded by Centrifuger, so the pattern must be quoted to keep the shell from expanding it first.
