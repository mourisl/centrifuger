---
title: Building an index
description: Turn reference FASTA files and an NCBI-style taxonomy into a Centrifuger index with centrifuger-build.
---

`centrifuger-build` compresses a collection of reference sequences into the FM-index that `centrifuger` searches. You only need to do this when no [pre-built index](/guides/prebuilt-indexes/) covers your reference set — for a custom genome collection, an in-house panel, or a database release that has not been published yet.

## What you need

Building an index requires three kinds of input:

| Input | Option | What it is |
|-------|--------|-----------|
| Reference sequences | `-r` or `-l` | FASTA files (optionally gzipped) holding the genomes or proteins |
| Taxonomy tree | `--taxonomy-tree` | NCBI `nodes.dmp` |
| Taxonomy names | `--name-table` | NCBI `names.dmp` |
| Sequence-to-taxon map | `--conversion-table` | Two columns: sequence ID, taxonomy ID |

Sequences can be given one file at a time with repeated `-r` options, or — much more practical for a large database — as a list file with `-l`, one path per line. The list file may carry the taxonomy ID in a second column, in which case `--conversion-table` is not needed.

:::note Different from Centrifuge
Centrifuge requires the reference files to be concatenated into a single FASTA. Centrifuger takes a *list* of files instead, so there is no giant intermediate file to create or store.
:::

## A minimal build

```bash
centrifuger-build -r ref.fa \
  --taxonomy-tree nodes.dmp \
  --name-table names.dmp \
  --conversion-table ref_seqid.map \
  -o cfr_ref_idx
```

This writes `cfr_ref_idx.1.cfr` through `cfr_ref_idx.4.cfr`. Those four files *are* the index — everything else used during the build can be deleted afterwards.

## Building a RefSeq index from scratch

`centrifuger-download` fetches both the taxonomy and the sequences from NCBI. The two commands below put the NCBI taxonomy in `taxonomy/` and all complete archaeal, bacterial and viral genomes in `library/`, capturing the sequence-to-taxon mapping as it goes:

```bash
centrifuger-download -o taxonomy taxonomy
centrifuger-download -o library -d "archaea,bacteria,viral" refseq > seqid2taxid.map
```

### Adding host or contaminant genomes

Including the host genome lets Centrifuger absorb host reads rather than mis-assigning them. Append to the same map file:

```bash
# human: T2T-CHM13
centrifuger-download -o library -d "vertebrate_mammalian" -t 9606 refseq >> seqid2taxid.map

# human: hg38 reference genome
centrifuger-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 \
  -c 'reference genome' refseq >> seqid2taxid.map

# mouse
centrifuger-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 10090 \
  -c 'reference genome' refseq >> seqid2taxid.map
```

`-t` selects a taxonomy ID (9606 is human, 10090 is mouse), `-a` restricts the assembly level, and `-c` restricts the RefSeq category.

### Collecting the files and building

Put the downloaded FASTA files into a list and hand it to `centrifuger-build`:

```bash
find library -type f -name "*.fna.gz" > file.list

centrifuger-build -t 4 \
  --conversion-table seqid2taxid.map \
  --taxonomy-tree taxonomy/nodes.dmp \
  --name-table taxonomy/names.dmp \
  -l file.list \
  -o refseq_abv \
  --build-mem 240G
```

:::tip If you used dustmasker
When `centrifuger-download` was run with low-complexity masking, the masked files are named `*_dustmasked.fna.gz`. Point `find` at that pattern instead.
:::

Once the build finishes, everything except `refseq_abv.[1234].cfr` can be removed.

## Controlling memory during the build

Index construction, not classification, is the memory-hungry step. Two parameters govern it:

- `--bmax INT` — block size for blockwise suffix-array sorting (default 16777216).
- `--dcv INT` — difference cover period (default 4096).

The defaults are tuned for modest databases and become inefficient for large ones. Rather than tuning them by hand, tell Centrifuger how much memory it may use and let it infer both:

```bash
centrifuger-build ... --build-mem 240G
```

`--build-mem` accepts `T`, `G`, `M` and `K` suffixes. Give it a realistic estimate of the memory actually available to the job — on a shared cluster, that is your allocation, not the machine's total.

:::caution Large builds take a long time
Indexing a database the size of RefSeq or `core_nt` is an hours-to-days job. Request a generous wall-clock limit, and consider `--checkpoint`.
:::

### Resuming an interrupted build

`--checkpoint` writes resumable state to `<output_prefix>_checkpoint.[123]`:

```bash
centrifuger-build ... --checkpoint
```

If the job is killed — by a wall-clock limit, a node failure or an out-of-memory event — rerunning the same command picks up from the last checkpoint instead of starting over.

## Shaping the contents of the index

A few options change *what* goes into the index rather than how it is built:

- `--subset-tax INT` — keep only the input genomes that sit under the given taxonomy node. Useful for building, say, a bacteria-only index from a broader download.
- `--concat-tax-genome` — concatenate all genomes sharing a taxonomy ID and drop the sequence ID information. This shrinks the index and speeds up search, at the cost of no longer knowing *which* sequence of a taxon a read matched.
- `--ignore-uncategorized-genome` — skip genomes whose sequence ID or taxonomy ID is missing or uncategorised. By default every input genome is included.
- `--offrate INT` — sample the suffix array every 2^INT BWT characters (default 4). A larger value makes the index smaller but resolving positions slower.

## Protein databases

Passing `--protein` builds an index over protein sequences instead of genomes:

```bash
centrifuger-build --protein -l protein_file.list \
  --taxonomy-tree taxonomy/nodes.dmp \
  --name-table taxonomy/names.dmp \
  --conversion-table seqid2taxid.map \
  -o cfr_protein_idx
```

See [protein classification](/guides/protein/) for how such an index is searched, and for the published protein indexes you can download instead.

## Other reference sources

The index-building procedure follows [Centrifuge's](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building) closely, with the program names changed (`centrifuge-download` → `centrifuger-download`, and so on), so recipes written for Centrifuge usually translate directly.

The `indices/` folder in the Centrifuger repository contains information for creating indexes from other sources, including **SILVA** and **GTDB**.

## After the build

- Keep only the `*.cfr` files; the downloaded FASTA and intermediate files are no longer needed for classification.
- The index prefix — the path without `.1.cfr` — is what you pass to `centrifuger -x` and `centrifuger-quant -x`.
- Indexes are portable: build once, copy to wherever the analysis runs.
