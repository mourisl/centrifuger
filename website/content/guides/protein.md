---
title: Protein classification
description: Build and search protein indexes with Centrifuger for translated classification against RefSeq proteins or NCBI nr.
---

Alongside nucleotide databases, Centrifuger can index and search **protein** sequences. Protein-level search is the tool of choice when the organisms in a sample are too divergent from anything in a nucleotide database for DNA matches to survive — deep-branching or poorly sampled taxa, and viruses in particular, where amino-acid sequence is conserved long after nucleotide sequence has drifted apart.

## Download a protein index

Two protein indexes are published; see [pre-built indexes](/guides/prebuilt-indexes/) for the full table.

| Title | Description | Size/~Memory | Date |
|-------|-------------|--------------|------|
| `cfr_protein_pv` | RefSeq bacteria, archaea and virus proteins | 25G | 2025/08/25 |
| `cfr_nr` | NCBI nr | 181G | 2026/01/31 |

```bash
centrifuger-download cfr_protein_pv
```

`cfr_protein_pv` is notably modest at 25 GB — smaller than any of the comprehensive nucleotide indexes — which makes translated search practical on hardware that could not hold `core_nt`.

## Build your own protein index

Pass `--protein` to `centrifuger-build`. Everything else works as it does for genomes: reference files via `-r` or `-l`, an NCBI-style taxonomy, and a sequence-to-taxon mapping.

```bash
centrifuger-build --protein \
  -l protein_files.list \
  --taxonomy-tree taxonomy/nodes.dmp \
  --name-table taxonomy/names.dmp \
  --conversion-table protein_seqid2taxid.map \
  -o cfr_protein_idx \
  -t 8 --build-mem 120G
```

The `--conversion-table` maps each protein accession to the taxonomy ID of the organism it came from, exactly as the nucleotide case maps contig accessions.

:::note The reference is protein, the reads are not
`--protein` describes the *reference*. You still supply ordinary nucleotide reads at classification time.
:::

See [building an index](/guides/building-an-index/) for `--build-mem`, `--checkpoint` and the other build-time options, all of which apply here too.

## Classify against it

Classification and quantification use the same commands as for a nucleotide index — only `-x` changes:

```bash
centrifuger -x cfr_protein_pv \
  -1 sample_1.fq.gz -2 sample_2.fq.gz -t 16 > protein_classification.tsv

centrifuger-quant -x cfr_protein_pv \
  -c protein_classification.tsv > protein_report.tsv
```

The output columns are unchanged, described in [output formats](/reference/output-formats/). The `seqID` column now names a protein accession rather than a contig, and `hitLength` reflects the matched extent of the read.

:::caution Scores are not comparable across alphabets
A score from a protein index and a score from a nucleotide index are on different scales. Compare protein runs with protein runs, and choose `--min-score` thresholds separately for each.
:::

## Nucleotide or protein?

| | Nucleotide index | Protein index |
|---|---|---|
| Best at | Species- and strain-level resolution among well-represented organisms | Detecting divergent or poorly sampled organisms |
| Resolution | Down to strain, given a strain-resolved database | Coarser — conserved proteins are shared across related taxa |
| Coverage of the read | Whole read, including non-coding sequence | Coding sequence only |
| Typical index size | 41–301 GB | 25–181 GB |

The two are complementary rather than competing. A common pattern is to classify against a nucleotide database first, then take the reads that came back unclassified and put them through a protein index:

```bash
centrifuger -x cfr_hpv+gbsarscov2 -1 s_1.fq.gz -2 s_2.fq.gz \
  --un unclassified -t 16 > nt_classification.tsv

centrifuger -x cfr_protein_pv \
  -1 unclassified_1.fq.gz -2 unclassified_2.fq.gz \
  -t 16 > protein_classification.tsv
```

That keeps the fast, high-resolution nucleotide search in front and reserves translated search for the residue that needs it.
