---
title: Citation and support
description: How to cite Centrifuger, where the source lives, and how to get help.
---

## Citing Centrifuger

If Centrifuger contributes to work you publish, please cite the Genome Biology paper:

> Song, L., Langmead, B. **Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification.** *Genome Biology* 25, 106 (2024). doi:[10.1186/s13059-024-03244-4](https://doi.org/10.1186/s13059-024-03244-4)

The paper received the **Best Paper Award at RECOMB 2024**.

- [Genome Biology, full text](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03244-4)
- [PubMed entry](https://pubmed.ncbi.nlm.nih.gov/38664753/) (PMID 38664753)
- [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.11.15.567129)

### BibTeX

```bibtex
@article{song2024centrifuger,
  title   = {Centrifuger: lossless compression of microbial genomes for
             efficient and accurate metagenomic sequence classification},
  author  = {Song, Li and Langmead, Ben},
  journal = {Genome Biology},
  volume  = {25},
  number  = {1},
  pages   = {106},
  year    = {2024},
  doi     = {10.1186/s13059-024-03244-4}
}
```

Please also cite the source of the reference database you used — RefSeq, GTDB, NCBI `nt`/`nr` and the LLNL-curated `core_nt` (PMID 40111052) each have their own citation.

## Support

Questions, bug reports and feature requests belong on the issue tracker:

**[github.com/mourisl/centrifuger/issues](https://github.com/mourisl/centrifuger/issues)**

A report is much easier to act on when it includes:

- the Centrifuger version and how it was installed (Bioconda or `make` from source);
- the exact command you ran;
- which index you used, and whether it was downloaded or self-built;
- the first lines of the error output, and any message on standard error;
- for classification problems, a few representative rows of the output.

The [FAQ and troubleshooting](/reference/faq/) page covers the questions that come up most often; it is worth a look before filing.

## Source code and licence

Centrifuger is developed in the open at [github.com/mourisl/centrifuger](https://github.com/mourisl/centrifuger).

Copyright © 2023–present, Li Song. The licence terms are in the `LICENSE` file of the repository.

## Related projects

- [Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/) — the predecessor whose index-building conventions and default report format Centrifuger follows.
- [GTDB](https://gtdb.ecogenomic.org/) — the Genome Taxonomy Database, source of the `cfr_gtdb_r232` indexes.
- [Bioconda](https://anaconda.org/bioconda/centrifuger) — the packaged distribution.

## About this documentation

This site is generated from Markdown sources in the `website/content` directory. Every page carries an **Edit this page on GitHub** link at the bottom; corrections and additions are welcome through the same repository as the code.
