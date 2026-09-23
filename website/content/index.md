---
title: Centrifuger — taxonomic classification on a compressed reference database
description: Fast, accurate taxonomic classification of sequencing reads against a losslessly compressed microbial genome or protein database.
---

<section class="hero"><div class="hero-inner"><h1>Centrifuger</h1><p class="tagline">Fast, accurate taxonmic classification and quantification using a <strong>losslessly compressed</strong> FM-index of comprehensive genome or protein sequence databases.</p><div class="cta"><a class="btn primary" href="getting-started/introduction/">Get started</a><a class="btn ghost" href="guides/prebuilt-indexes/">Download an index</a><a class="btn ghost" href="https://github.com/mourisl/centrifuger" target="_blank" rel="noopener">View on GitHub</a></div></div></section>

<section class="home-section">

## What is Centrifuger?

Centrifuger is an efficient taxonomic classification method that compares sequencing reads against a microbial genome or protein database. It implements a novel lossless compression method, the **run-block compressed BWT**, together with other strategies for compacting the Ferragina–Manzini (FM) index. The result is a database representation that is small enough to hold in the memory of an ordinary server, yet loses none of the underlying sequence.

For example, Centrifuger can classify reads against the 2023 RefSeq prokaryotic genomes — about 140 billion nucleotides — using roughly **43 GB of memory**. Because the compressed index is still directly searchable, classification does not have to be slow: a typical sequencing sample is processed within an hour.

<div class="stats"><div class="stat"><div class="value">140 Gbp</div><div class="label">RefSeq prokaryotic genomes indexed</div></div><div class="stat"><div class="value">43 GB</div><div class="label">Memory to classify against that database</div></div><div class="stat"><div class="value">~2&times;</div><div class="label">Smaller than other FM-index classifiers</div></div><div class="stat"><div class="value">Lossless</div><div class="label">No k-mer or minimizer subsampling</div></div></div>

</section>

<section class="home-section">

## Why Centrifuger

<div class="cards"><div class="card"><h3>Lossless compression</h3><p>The run-block compressed BWT and a hybrid run-length compressed BWT exploit the intermediate repetitiveness of microbial genome collections, reaching sublinear space without discarding sequence.</p></div><div class="card"><h3>Better species-level calls</h3><p>Because nothing is subsampled and match length is not capped by a fixed <em>k</em>, Centrifuger improves both sensitivity and precision at the species and genus level.</p></div><div class="card"><h3>Indexes ready to download</h3><p>Pre-built indexes for RefSeq, GTDB r232, NCBI core nt and NCBI nr are published, and <code>centrifuger-download</code> fetches them for you.</p></div><div class="card"><h3>Short reads, long reads, single cells</h3><p>Paired-end, single-end, interleaved and long-read input all work, and barcode/UMI parsing supports 10x Genomics and combinatorial-barcoding protocols such as SHARE-seq.</p></div><div class="card"><h3>Nucleotide and protein search</h3><p>Build a protein index with <code>--protein</code> to classify reads through translated search against RefSeq proteins or NCBI nr.</p></div><div class="card"><h3>Reports other tools understand</h3><p><code>centrifuger-quant</code> writes abundance profiles in the Centrifuge, MetaPhlAn, CAMI or Kraken-report format.</p></div></div>

</section>

<section class="quickstart"><div class="home-section">

## Quick start

Install from Bioconda, download a pre-built index, and classify a paired-end sample:

```bash
conda install -c conda-forge -c bioconda centrifuger

centrifuger-download cfr_hpv+gbsarscov2

centrifuger -x cfr_hpv+gbsarscov2 -1 sample_1.fq.gz -2 sample_2.fq.gz -t 8 \
  > classification.tsv

centrifuger-quant -x cfr_hpv+gbsarscov2 -c classification.tsv > report.tsv
```

The full walkthrough, including building your own index, is in the [quick start guide](/getting-started/quick-start/).

</div></section>

<section class="home-section">

## Citation

<div class="citation"><p>Song, L., Langmead, B. <strong>Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification.</strong> <em>Genome Biology</em> 25, 106 (2024).</p><p><a href="https://doi.org/10.1186/s13059-024-03244-4" target="_blank" rel="noopener">doi:10.1186/s13059-024-03244-4</a> &middot; <a href="https://pubmed.ncbi.nlm.nih.gov/38664753/" target="_blank" rel="noopener">PubMed</a></p><p><span class="award">Best Paper Award at RECOMB 2024</span></p></div>

Centrifuger is copyright &copy; 2023–present, Li Song. Questions and bug reports are welcome on the [GitHub issue tracker](https://github.com/mourisl/centrifuger/issues).

</section>
