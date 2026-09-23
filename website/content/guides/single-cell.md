---
title: Single-cell and barcoded data
description: Classify 10x Genomics and other barcoded libraries, including UMI extraction, whitelists, barcode translation and reading barcodes from BAM headers.
---

Centrifuger can carry cell barcodes and UMIs through classification, so that per-read assignments can be grouped by cell afterwards. This makes it usable for microbial detection in single-cell libraries — 10x Genomics, combinatorial barcoding schemes such as SHARE-seq, and anything else whose barcode position you can describe.

## The pieces

| Option | Purpose |
|--------|---------|
| `--barcode STR` | File containing the cell barcode sequence |
| `--UMI STR` | File containing the UMI sequence |
| `--read-format STR` | Where in those files read, barcode and UMI actually live |
| `--barcode-whitelist STR` | List of valid barcodes, for error correction |
| `--barcode-translate STR` | Map observed barcodes onto another set of identifiers |

If `--barcode` or `--UMI` is omitted while `--read-format` asks for one, Centrifuger extracts the pattern from read 1.

:::note Whitelists need an explicit barcode file
`--barcode-whitelist` requires `--barcode`. It has nothing to look up otherwise.
:::

## The --read-format specification

`--read-format` is a comma-separated list of fields. Each field is colon-separated:

```text
[r1|r2|bc|um]:start:end:strand
```

- The first part names what is being extracted: read 1, read 2, the barcode, or the UMI.
- `start` and `end` are **0-based and inclusive**. `-1` means "to the end of the read".
- `strand` is `+` or `-`. With `-`, the extracted sequence is reverse-complemented. It may be omitted when it is `+`, and it is ignored for `r1` and `r2`.

Multiple fields of the same kind describe non-consecutive segments, which are concatenated in order — `bc:0:15,bc:32:-1` takes bases 0–15 and 32 onward as one barcode.

### Example: barcode in the first 16 bp of read 1

The barcode occupies the first 16 bases of read 1, and the biological sequence is everything after it:

```bash
centrifuger -x cfr_idx \
  -1 read1.fq.gz -2 read2.fq.gz \
  --barcode read1.fq.gz \
  --read-format bc:0:15,r1:16:-1
```

Note that `read1.fq.gz` is passed twice: once as a read file and once as the barcode file. `--read-format` then splits it.

## 10x Genomics libraries

In a typical 10x 3' library, R1 carries a 16 bp cell barcode followed by a 12 bp UMI, and R2 carries the cDNA. Classify R2 while pulling barcode and UMI out of R1:

```bash
centrifuger -x cfr_idx \
  -u "path_to_10x_fastqs/*_R2_*.fastq.gz" \
  --barcode "path_to_10x_fastqs/*_R1_*.fastq.gz" \
  --UMI "path_to_10x_fastqs/*_R1_*.fastq.gz" \
  --read-format bc:0:15,um:16:-1 \
  --barcode-whitelist cellranger_folder/cellranger-cs/VERSION/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
```

:::caution Two things to get right
The exact offsets depend on your 10x chemistry — check the kit before copying the numbers above. And the quotes around the wildcard paths are **necessary**, so that Centrifuger expands the pattern rather than the shell.
:::

## Reading barcodes from FASTQ header comments

When reads have already been through an aligner, the barcode often lives in a tag in the read header rather than in the sequence. Centrifuger can parse it from there using an extended form:

```text
[bc|um]:hd:field:start:end:strand
```

`hd` is a keyword meaning "search the header comment". `field` is either

- a **number** — the 0-based index of the whitespace-separated field in the comment, not counting the read ID; or
- a **string** — a prefix to search for, with extraction starting from the match.

For a header like:

```text
@r1 CR:Z:NNNN CB:Z:ACGT UR:Z:NNNN
```

both of these extract the corrected barcode `ACGT`:

```text
bc:hd:1:5:-1      # field 1 of the comment, from offset 5 to the end
bc:hd:CB:5:-1     # the field starting with "CB", from offset 5 to the end
```

The offset of 5 skips the `CB:Z:` tag prefix.

### From a Cell Ranger BAM

This is what makes it practical to look for microbial sequence among the reads that did not map to the host. Single-end:

```bash
samtools fastq -T CB -f 4 alignment.bam \
  | centrifuger -x cfr_idx --read-format bc:hd:CB:5:-1 -u - > output.tsv
```

Paired-end, keeping both barcode and UMI and saving whatever stays unclassified:

```bash
samtools view -f 0xC alignment.bam \
  | samtools sort -n - \
  | samtools fastq -T CB,UB - \
  | centrifuger -x cfr_idx \
      --read-format bc:hd:CB:5:-1,um:hd:UB:5:-1 \
      --un unclassified -i - > output.tsv
```

`-f 4` keeps unmapped reads; `-f 0xC` keeps fragments where both mates are unmapped. `samtools fastq -T` copies the named tags into the FASTQ header, which is what `bc:hd:...` then reads. Sorting by name (`samtools sort -n`) is what allows the two mates to be emitted as an interleaved stream for `-i -`.

## Barcode translation and combinatorial barcoding

`--barcode-translate` maps observed barcodes onto another set of identifiers. The translation file is a two-column TSV or CSV:

```text
translated_barcode    original_barcode
```

— the translated identifier first, the observed sequence second.

This also covers **combinatorial barcoding** such as SHARE-seq, where a cell is identified by several barcode segments rather than one sequence. Centrifuger translates each segment listed in the second column to the identifier in the first, and joins the resulting identifiers with `-` in the output. Combine it with a multi-segment `--read-format`:

```bash
centrifuger -x cfr_idx -u reads.fq.gz \
  --barcode barcodes.fq.gz \
  --read-format bc:0:7,bc:38:45,bc:76:83 \
  --barcode-translate round_barcodes.tsv
```

## Putting the output to use

Barcode and UMI columns ride along in the classification file, so grouping by cell afterwards is ordinary table work: count distinct UMIs per (barcode, taxID) pair to get a per-cell microbial profile, and filter to the barcodes your single-cell pipeline called as real cells.

Passing `--barcode-whitelist` is worth the effort here — correcting sequencing errors in the barcode keeps reads from a single cell from being scattered across dozens of phantom barcodes.

## Next steps

- [Classifying reads](/guides/classification/) — the rest of the classification options.
- [Command-line interface](/reference/cli/#centrifuger) — the full option list.
