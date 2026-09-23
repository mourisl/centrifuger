---
title: Pre-built indexes
description: Published Centrifuger indexes for RefSeq, GTDB, NCBI core nt and NCBI nr, and how to download them.
---

Building a comprehensive index is an hours-to-days job on a large machine. For most projects there is no need: the Centrifuger authors publish ready-made indexes for the databases people ask for most often. Download one, point `-x` at it, and start classifying.

## Available indexes

The **Size/~Memory** column is both the download size and a good estimate of the RAM needed to classify against the index.

| Title | Description | Link | Size/~Memory | Date |
|-------|-------------|------|--------------|------|
| **Genome databases** | | | | |
| `cfr_hpv+gbsarscov2` | RefSeq human, bacteria, archaea, virus + SARS-CoV-2 variants from GenBank | [Zenodo](https://zenodo.org/records/10023239) | 41G | 2023/10/01 |
| `cfr_gtdb_r232` | [GTDB r232](https://gtdb.ecogenomic.org/) | [Dropbox](https://www.dropbox.com/scl/fo/wilot7fzf8bzyfsco81wi/AFzwZhSgE9N1r7uai8DeUec?rlkey=yazr6e4b71pxbz4o4bu7dqft5&st=ekr030j8&dl=0) | 230G | 2026/05/15 |
| `cfr_gtdb_r232+refseq_hvfpc` | GTDB r232 + RefSeq human, virus, fungi, protozoa and contaminants (UniVec, EmVec) | [Dropbox](https://www.dropbox.com/scl/fo/6sf1k7iwnckdca2opl1ir/AK8X0gPn90TF7Ox0XxK8kjI?rlkey=y7787xqylcfssvwjn1w66rmqz&st=ck0sa631&dl=0) | 232G | 2026/05/15 |
| `cfr_core_nt` | NCBI core nt | [Dropbox](https://www.dropbox.com/scl/fo/f1mbf7nf893pisoruanb4/AHS06LaJr9EN0Pg7hbifWn8?rlkey=7fgtj6pi53l2iwrjw1k6xq8o8&st=yn57lnkh&dl=0) | 212G | 2025/06/11 |
| `cfr_llnl_core_nt_202603` | LLNL-curated NCBI core_nt (PMID:40111052) | [Dropbox](https://www.dropbox.com/scl/fo/zkjoh4luk5e3kvzs9hvir/APWYWaOHu-7RVxELTi41_rI?rlkey=soavq4el3od1nch7op99n5mz5&st=o1ckosy7&dl=0) | 242G | 2026/03/01 |
| `cfr_llnl_core_nt_wseqid_202512` | LLNL-curated NCBI core_nt with sequence ID information but less synchronised taxdmp | [Dropbox](https://www.dropbox.com/scl/fo/zkjoh4luk5e3kvzs9hvir/APWYWaOHu-7RVxELTi41_rI?rlkey=soavq4el3od1nch7op99n5mz5&st=o1ckosy7&dl=0) | 301G | 2025/12/01 |
| **Protein databases** | | | | |
| `cfr_protein_pv` | RefSeq bacteria, archaea and virus proteins | [Zenodo](https://zenodo.org/records/22663514) | 25G | 2025/08/25 |
| `cfr_nr` | NCBI nr | [Dropbox](https://www.dropbox.com/scl/fo/nfnm3nehfmx3or3anrvnk/ABqABdp5-LZP_AtD9Y0zXWE?rlkey=qt87y5966jvv60s7gcxjsha2q&st=m71m0buz&dl=0) | 181G | 2026/01/31 |

Older indexes remain available in [this Dropbox folder](https://www.dropbox.com/scl/fo/08horwj8mdzarlk2ocyky/AJIUqBg4ZU4qXdaTBnl64xM?rlkey=y7vk78c3o1pd2fq20f258vuyf&st=57xyuyjl&dl=0).

:::note This table mirrors the repository
The canonical list lives in the [Centrifuger README](https://github.com/mourisl/centrifuger#build-index) and gains entries as new database releases are indexed. Check there if you need something newer than what is shown above.
:::

## Downloading with centrifuger-download

The simplest route is to let Centrifuger fetch the index by title:

```bash
centrifuger-download cfr_hpv+gbsarscov2
```

Substitute any title from the first column of the table.

## Downloading by hand

For the Dropbox-hosted indexes you may prefer to pull the files yourself — for instance on a cluster node where you want to control where the download lands, or to resume a partial transfer. Open the Dropbox folder in a browser, right-click each file and choose **Copy link**, then fetch it with `wget`:

```bash
wget -O cfr_core_nt.1.cfr "https://www.dropbox.com/...&dl=1"
```

:::tip Force a direct download
Dropbox share links end in `dl=0`, which serves an HTML preview page. Change it to `dl=1` so `wget` or `curl` receives the file itself.
:::

Download every `*.cfr` file belonging to the index and keep them together in one directory with a common prefix.

## Choosing an index

| If you want to… | Consider |
|-----------------|----------|
| Classify human-associated microbiome samples on a mid-sized server | `cfr_hpv+gbsarscov2` (41G) |
| Use a standardised, rank-normalised bacterial and archaeal taxonomy | `cfr_gtdb_r232` |
| Cover prokaryotes *and* host, fungi, protozoa and vector contaminants | `cfr_gtdb_r232+refseq_hvfpc` |
| Search the broadest nucleotide collection, including eukaryotes | `cfr_core_nt` or an LLNL-curated variant |
| Detect divergent organisms through translated search | `cfr_protein_pv` or `cfr_nr` |

Two practical considerations tend to decide it:

- **Memory.** The index must fit in RAM. A 230 GB index needs a machine with well over 230 GB — the 41 GB RefSeq index runs comfortably where the large ones cannot run at all.
- **Host reads.** If your sample contains host DNA, prefer an index that includes the host genome. Reads that have nowhere correct to go are the main source of spurious microbial calls.

## Using a downloaded index

Pass the prefix — the path without the `.1.cfr` suffix — to `-x`:

```bash
centrifuger -x /data/indexes/cfr_hpv+gbsarscov2 \
  -1 sample_1.fq.gz -2 sample_2.fq.gz -t 16 > classification.tsv

centrifuger-quant -x /data/indexes/cfr_hpv+gbsarscov2 \
  -c classification.tsv > report.tsv
```

:::caution Fast storage helps the first load
The index is read into memory at startup. On a shared filesystem, loading a 200 GB index can take longer than the classification itself; a local scratch disk makes a noticeable difference when you run many samples.
:::

If none of these databases fits your project, [build your own](/guides/building-an-index/).
