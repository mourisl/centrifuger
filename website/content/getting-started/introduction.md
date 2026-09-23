---
title: Introduction
description: What Centrifuger is, how it compresses a reference database, and when to reach for it.
---

Centrifuger is an efficient taxonomic classification method that compares sequencing reads against a microbial genome or protein database. Given a set of reads and a reference database such as RefSeq, GTDB or NCBI `nr`, Centrifuger reports, for every read, which reference sequence it came from and the taxonomic identifier of that sequence. A companion program turns those per-read assignments into an abundance profile for the sample.

What distinguishes Centrifuger from other classifiers is how the reference database is stored. It implements a novel lossless compression method — the **run-block compressed BWT** — and combines it with other strategies for compacting the Ferragina–Manzini (FM) index. The database stays complete: no *k*-mers are dropped, no minimizers are subsampled, and no sequence is discarded. It simply takes much less room.

## Why compression matters

Comprehensive reference databases have grown faster than the memory of the machines used to search them. The 2023 RefSeq prokaryotic collection is roughly 140 billion nucleotides; NCBI `core_nt` and `nr` are larger still. Classifiers generally cope in one of two ways: they throw information away — subsampling *k*-mers, restricting matches to a fixed length — or they demand a very large machine.

Centrifuger takes the third route. Microbial genome databases sit at an *intermediate* level of repetitiveness: far more repetitive than a single genome, far less repetitive than a collection of thousands of strains of one species. Centrifuger introduces two compact data structures tuned for exactly that regime:

- **Run-block compressed BWT (RBBWT)** — a run-block compressed sequence representation that achieves sublinear storage for the BWT without sacrificing much time efficiency.
- **Hybrid run-length compressed BWT** — a complementary representation used where run-length coding pays off.

Together with other FM-index compaction strategies, this halves the memory footprint compared with other FM-index-based approaches, letting the 2023 RefSeq prokaryotic database be classified against in about 43 GB of memory.

## Why lossless matters

Compression is not only an engineering convenience here — it is what makes the accuracy possible.

Because the representation is lossless and the index is a true FM-index, Centrifuger can extend a match for as long as the read and the reference agree. There is no fixed *k* that truncates a match, and no sampled *k*-mer set that can miss one. That unconstrained match length, together with the complete database, is what gives Centrifuger significantly better sensitivity **and** precision for species- and genus-level classification.

:::tip What you gain in practice
Longer, exactly resolved matches make it easier to tell apart closely related genomes, which is precisely where species-level metagenomic classification is hardest.
:::

## What is in the distribution

Centrifuger ships as a small set of command-line programs:

| Program | Purpose |
|---------|---------|
| `centrifuger-build` | Build a searchable index from reference FASTA files and an NCBI-style taxonomy |
| `centrifuger` | Classify reads against an index; writes a per-read TSV to standard output |
| `centrifuger-quant` | Turn a classification file into an abundance profile |
| `centrifuger-download` | Download reference sequences and taxonomy from NCBI, or fetch a pre-built index |

The usual workflow is a straight line through them:

```text
reference FASTA + taxonomy  ──centrifuger-build──▶  index (*.cfr)
                                                       │
                       reads (FASTQ)  ──centrifuger────┘──▶  classification.tsv
                                                                    │
                                        centrifuger-quant ──────────┘──▶  report.tsv
```

If a published index already covers your reference set, you can skip the first step entirely — see [pre-built indexes](/guides/prebuilt-indexes/).

## What Centrifuger handles

- **Read layouts** — paired-end (`-1`/`-2`), single-end (`-u`), interleaved (`-i`), and many samples at once through a sample sheet. The option follows the layout of the data, not its read length.
- **Read lengths** — short-read and long-read data, with no separate mode for either.
- **Assay types** — bulk sequencing as well as barcoded single-cell data, including 10x Genomics and combinatorial barcoding such as SHARE-seq.
- **Sequence alphabets** — nucleotide databases, and protein databases through `centrifuger-build --protein`.
- **Report formats** — the Centrifuge report format plus MetaPhlAn, CAMI and Kraken-report output from `centrifuger-quant`.

## Relationship to Centrifuge

Centrifuger follows the design and conventions of [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/), from the same group, closely enough that the index-building procedure and the default report format will look familiar; the programs are renamed (`centrifuge-download` becomes `centrifuger-download`, and so on). The important differences are the compressed index, the unconstrained match length, and a few conveniences such as file-list input instead of concatenated FASTA. If you are coming from Centrifuge, [building an index](/guides/building-an-index/) points out where the two diverge.

## Where to go next

- [Installation](/getting-started/installation/) — Bioconda, or `make` from source.
- [Quick start](/getting-started/quick-start/) — classify a sample end to end with the bundled example data.
- [Command-line interface](/reference/cli/) — every option of every program in one place.
