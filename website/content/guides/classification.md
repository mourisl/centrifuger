---
title: Classifying reads
description: Run centrifuger on paired-end, single-end, interleaved, long-read and multi-sample input, and tune how hits are reported.
---

`centrifuger` is the classification program. It loads an index, streams reads past it, and writes one row per read assignment to standard output.

```bash
centrifuger -x cfr_idx -1 read_1.fq.gz -2 read_2.fq.gz -t 8 > output.tsv
```

Output goes to stdout, so redirect it to a file or pipe it onward.

## Specifying reads

Exactly one input form is required, alongside `-x`.

| Form | Option | Example |
|------|--------|---------|
| Paired-end | `-1 FILE -2 FILE` | `-1 s_1.fq.gz -2 s_2.fq.gz` |
| Single-end | `-u FILE` | `-u reads.fq.gz` |
| Interleaved paired-end | `-i FILE` | `-i interleaved.fq.gz` |
| Many samples at once | `--sample-sheet FILE` | `--sample-sheet samples.tsv` |

Input may be gzipped or uncompressed; Centrifuger detects which.

### Reading from a pipe

Use `-` as the filename to read from standard input. This is how you classify reads that are produced by another program rather than stored on disk — for example the unmapped reads of an alignment:

```bash
samtools fastq -f 4 alignment.bam | centrifuger -x cfr_idx -u - > output.tsv
```

### Wildcards

`-1`, `-2` and `-u` accept shell wildcards, which is convenient for sequencer output split across lanes. **Quote the pattern** so that Centrifuger, and not the shell, expands it:

```bash
centrifuger -x cfr_idx -u "run/*_R2_*.fastq.gz" > output.tsv
```

:::caution The quotes are required
Without them the shell expands the pattern into many arguments and Centrifuger sees only the first one.
:::

### Sample sheets

To classify many samples in a single invocation, list them in a sample sheet. Each row has five fields:

```text
read1   read2   barcode   UMI   output
```

Use a dot (`.`) for any field that does not apply. For plain paired-end samples with no barcode or UMI:

```text
sampleA_1.fq.gz  sampleA_2.fq.gz  .  .  sampleA.tsv
sampleB_1.fq.gz  sampleB_2.fq.gz  .  .  sampleB.tsv
```

```bash
centrifuger -x cfr_idx --sample-sheet samples.tsv -t 16
```

Because each row names its own output file, results are written per sample rather than to stdout. The index is loaded once and reused across every sample, which is the main reason to prefer a sample sheet over a shell loop when the index is large.

## Threads

`-t INT` sets the number of classification threads (default 1). Classification parallelises well, so set this to the cores available to your job:

```bash
centrifuger -x cfr_idx -1 s_1.fq.gz -2 s_2.fq.gz -t 16 > output.tsv
```

Threads share one copy of the index, so raising `-t` does not raise memory use appreciably.

## How many assignments per read

`-k INT` controls how many distinct, primary assignments are reported for each read or read pair. The default is 1.

```bash
centrifuger -x cfr_idx -1 s_1.fq.gz -2 s_2.fq.gz -k 5 > output.tsv
```

A read that matches several genomes equally well — common for conserved regions and for closely related strains — produces up to `k` rows. This matters downstream: [`centrifuger-quant`](/guides/quantification/) distributes multi-classified reads across the taxa they hit, so a larger `-k` gives a more ambiguous but more specific picture, and can improve the abundance estimate.

The related `--hitk-factor INT` caps the work spent resolving each hit, at `k × factor` entries (default 40; `0` removes the restriction). Raise it if you use a large `-k` on a highly redundant database and suspect assignments are being truncated.

## Tuning what counts as a hit

- `--min-hitlen INT` — minimum length of a partial hit. The default is chosen automatically from the data. Raising it makes classification stricter and reduces spurious short matches; lowering it recovers hits from very short or heavily degraded reads.
- `--consider-secondary INT,FLOAT` — accept a secondary hit when its hit length is at least `INT` and its score is at least `FLOAT` times the best score. The default is `2000,0.995`, i.e. only long, near-equal alternatives.
- `--no-dust` — disable DUST masking of low-complexity regions in reads. Masking is on by default and prevents homopolymers and simple repeats from generating misleading matches; turn it off only if you specifically need low-complexity sequence considered.
- `--merge-readpair` — merge overlapping paired-end reads and trim adapters before classification. Worth enabling for short inserts, where the merged fragment gives a longer, more specific match than either mate alone.

## Saving classified and unclassified reads

Two options write reads back out, using the given prefix and the input's layout (`<prefix>_1.fq.gz` / `<prefix>_2.fq.gz` for paired input):

```bash
centrifuger -x cfr_idx -1 s_1.fq.gz -2 s_2.fq.gz \
  --un unclassified --cl classified > output.tsv
```

- `--un STR` — reads that received no assignment. Useful for assessing database coverage, or for passing the remainder to assembly or a broader database.
- `--cl STR` — reads that were assigned. Useful for host depletion when the host genome is in the index: classify, keep the unclassified reads, and continue.

## Long reads

Long reads need no special flag. Pick the input option by layout, exactly as for short reads — long-read data is normally single-end, so it goes through `-u`:

```bash
centrifuger -x cfr_idx -u nanopore.fq.gz -t 16 > output.tsv
```

Because match length is not capped by a fixed *k*, long reads often produce very long exact matches and correspondingly high scores. If a run produces many equally good assignments across related genomes, raise `-k` so that the ambiguity is visible in the output rather than resolved arbitrarily.

## Worked examples

```bash
# single-end data
centrifuger -x cfr_idx -u read_single_end.fq > output.tsv

# paired-end data
centrifuger -x cfr_idx -1 read_1.fq -2 read_2.fq > output.tsv

# 10x Genomics single-end fastq files, with barcode and UMI
centrifuger -x cfr_idx -u "*_R2_*.fastq.gz" \
  --barcode "*_R1_*.fastq.gz" --UMI "*_R1_*.fastq.gz" \
  --read-format bc:0:15,um:16:-1 \
  --barcode-whitelist 3M-february-2018.txt.gz > output.tsv

# unmapped reads from a single-end Cell Ranger BAM
samtools fastq -T CB -f 4 alignment.bam \
  | centrifuger -x cfr_idx --read-format bc:hd:CB:5:-1 -u - > output.tsv

# unmapped fragments from a paired-end Cell Ranger BAM, keeping barcode and UMI
samtools view -f 0xC alignment.bam \
  | samtools sort -n - \
  | samtools fastq -T CB,UB - \
  | centrifuger -x cfr_idx \
      --read-format bc:hd:CB:5:-1,um:hd:UB:5:-1 \
      --un unclassified -i - > output.tsv
```

The barcode-aware invocations are explained in [single-cell and barcoded data](/guides/single-cell/).

## Next steps

- [Output formats](/reference/output-formats/#classification-output) — what each column of `output.tsv` means.
- [Abundance quantification](/guides/quantification/) — turning assignments into a profile.
- [Command-line interface](/reference/cli/#centrifuger) — the complete option list.
