---
title: FAQ and troubleshooting
description: Memory requirements, choosing a database, interpreting scores, and what to do when a run does not behave.
---

## Installation and setup

### How much memory do I need?

Enough to hold the index. The **Size/~Memory** column of the [pre-built index table](/guides/prebuilt-indexes/) is a good estimate for each published database:

- `cfr_hpv+gbsarscov2` (RefSeq human, bacteria, archaea, virus + SARS-CoV-2) — about 41 GB. Centrifuger can classify against the 2023 RefSeq prokaryotic genomes, roughly 140 billion nucleotides, in about 43 GB of memory.
- `cfr_protein_pv` (RefSeq prokaryotic and viral proteins) — about 25 GB, the smallest comprehensive option.
- `core_nt` and GTDB-scale indexes — 200–300 GB.

Threads share a single copy of the index, so raising `-t` barely changes memory use. Building an index is a different matter: see [why does my build run out of memory?](/reference/faq/#why-does-my-build-run-out-of-memory)

### Do I have to build an index?

Usually not. Published indexes cover RefSeq, GTDB r232, NCBI `core_nt` and NCBI `nr`, and `centrifuger-download <title>` fetches them. Build your own when you need a custom genome collection, an in-house panel, or a database release that has not been indexed yet.

### Can I reuse a Centrifuge index?

No. Centrifuger stores the database in its own run-block compressed FM-index and reads `*.cfr` files. The *procedure* for building one closely follows Centrifuge's, and `centrifuge-download` recipes translate directly to `centrifuger-download`, but the index itself must be built or downloaded for Centrifuger.

### Does an upgrade invalidate my index?

Not normally — indexes are portable across machines and survive routine upgrades. If a release ever changes the index format, the GitHub release notes will say so.

## Running classifications

### How do I choose -k?

`-k` is how many distinct, primary assignments are reported per read, default 1.

- Keep `-k 1` when you want a single best call per read and a compact output file.
- Raise it when reads are shared between closely related genomes. `centrifuger-quant` spreads a multi-classified read across the taxa it hit, so a larger `-k` yields an ambiguous but more specific classification, and can improve the abundance estimate.

If species-level proportions matter, classify once with `-k 1` and once with `-k 5` and compare the profiles.

### What is a good score cutoff?

There is no universal value — the score is a weighted sum of hits, so it scales with read length and with how much of the read matched. Rather than guessing, run without filters, inspect the distribution of column 4 in the classification file, and choose a threshold that removes the low-score tail:

```bash
tail -n +2 classification.tsv | cut -f4 | sort -n | uniq -c | tail -30
```

Then apply it at quantification time with `--min-score`, which does not require re-running the expensive classification step. `--min-length` does the same for classified length.

### My reads are long. Do I need different options?

No special flag. Read length does not change how you invoke Centrifuger — you choose the input option by the *layout* of the data, so long reads go through `-u` because they are single-end, not because they are long.

What does change is the result. Because match length is not capped by a fixed *k*, long reads tend to produce long matches and correspondingly high scores. If a run reports many equally good assignments across related genomes, raise `-k` so the ambiguity appears in the output rather than being resolved arbitrarily.

### Should I turn off DUST masking?

Rarely. Low-complexity masking is on by default and stops homopolymers and simple repeats from generating misleading matches. `--no-dust` disables it, which is appropriate only when low-complexity sequence is itself the object of study.

### Why do my wildcards only match one file?

The shell expanded them first. Centrifuger does its own wildcard expansion for `-1`, `-2` and `-u`, so the pattern must be quoted:

```bash
centrifuger -x cfr_idx -u "run/*_R2_*.fastq.gz" > output.tsv
```

### How do I classify many samples efficiently?

Use `--sample-sheet` rather than a shell loop. Each row is `read1 read2 barcode UMI output`, with a dot (`.`) for fields that do not apply, and the index is loaded once and reused for all samples. With a 200 GB index, that alone can dominate the total runtime of a batch.

## Building indexes

### Why does my build run out of memory?

The defaults for `--bmax` (16777216) and `--dcv` (4096) are tuned for modest databases and become inefficient for large ones. Do not tune them by hand — tell Centrifuger the memory budget and let it infer both:

```bash
centrifuger-build ... --build-mem 240G
```

Give a realistic figure: on a shared cluster that is your job's allocation, not the machine's total memory.

### My build was killed halfway through. Do I start over?

Not if you passed `--checkpoint`. It writes `<prefix>_checkpoint.[123]`, and rerunning the same command resumes from the last checkpoint. For a database of RefSeq or `core_nt` scale — an hours-to-days build — it is worth adding by default.

### Can I make the index smaller?

Several options trade something for size:

- `--concat-tax-genome` concatenates genomes sharing a taxID and drops seqID information. Smaller and faster to search, but you no longer learn *which* sequence of a taxon a read matched.
- `--subset-tax INT` keeps only genomes under one taxonomy node.
- `--ignore-uncategorized-genome` skips genomes with a missing or uncategorised seqID or taxID.
- `--offrate INT` samples the suffix array more sparsely: a smaller index, slower position resolution.

### Do I need to concatenate my FASTA files?

No — and this is a difference from Centrifuge. Centrifuger takes a list of files via `-l`, one path per row, with an optional taxonomy ID in a second column.

## Interpreting results

### Why is a taxon in my report that cannot be in my sample?

Some common causes, roughly in order of likelihood:

- **The host genome is not in the index.** Host reads with nowhere correct to go land on whatever is closest, which is the largest single source of spurious calls. Prefer an index that includes your host.
- **Low-complexity or adapter sequence.** Keep DUST masking on, and consider `--merge-readpair` for short inserts so adapters are trimmed.
- **A thin database.** If the true organism is absent, its reads are distributed among relatives.
- **Genuinely weak evidence.** Check `numUniqueReads` — a taxon with many reads but almost no unique ones is not well supported. Filter with `--min-score` and `--min-length`.

### numReads and numUniqueReads differ a lot. Which do I trust?

`numUniqueReads` is the conservative evidence: reads assigned unambiguously below that node. `numReads` includes multi-classified reads, distributed evenly among the taxa they hit. A wide gap means the taxon is difficult to separate from its neighbours, so treat its individual abundance with caution while the clade as a whole may be solid.

### Why does my report mix strains, species and genera?

Reads that cannot be resolved to a strain are attributed to the lowest rank that does explain them. A profile that mixes ranks is normal and is more honest than forcing every read down to strain level.

### How do I get a Kraken-style report?

```bash
centrifuger-quant -x cfr_idx -c classification.tsv --output-format 3 > kreport.txt
```

`--output-format` also emits MetaPhlAn (`1`) and CAMI (`2`) profiles from the same classification file.

### Almost nothing was classified. What now?

1. Confirm the index matches the sample — a prokaryote-only index will not classify a eukaryotic or viral sample well.
2. Look at the unclassified reads themselves: `--un unclassified` writes them out. If they are adapter, host or low-complexity sequence, that is your answer.
3. Try a broader database such as `cfr_core_nt`, or a protein index for [translated search](/guides/protein/) against divergent organisms.
4. Check that read quality and length are what you expect; `--min-hitlen` defaults to an automatic value, but severely truncated reads may still fall short.

## Getting help

If none of this covers your situation, open a [GitHub issue](https://github.com/mourisl/centrifuger/issues). Including your Centrifuger version, the exact command, the index you used, and the first lines of the error output makes a diagnosis much faster.
