---
title: Installation
description: Install Centrifuger from Bioconda or build it from source with make.
---

Centrifuger is distributed both as a Bioconda package and as source code that compiles with a single `make`. Either route gives you the same four executables: `centrifuger-build`, `centrifuger`, `centrifuger-quant` and `centrifuger-download`.

## Requirements

- A Linux or macOS system with a C++ compiler and `make`.
- [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads) — Centrifuger's only dependency, and present by default on essentially every Unix-like system.
- Enough RAM to hold the index you intend to use. This is the real constraint: see [how much memory do I need?](/reference/faq/#how-much-memory-do-i-need)

## Install with Conda (recommended)

Centrifuger is available from [Bioconda](https://anaconda.org/bioconda/centrifuger):

```bash
conda install -c conda-forge -c bioconda centrifuger
```

Installing into a dedicated environment keeps it isolated from the rest of your tooling:

```bash
conda create -n centrifuger -c conda-forge -c bioconda centrifuger
conda activate centrifuger
```

Mamba works as a drop-in replacement and resolves the environment considerably faster:

```bash
mamba install -c conda-forge -c bioconda centrifuger
```

## Build from source

Clone the [GitHub repository](https://github.com/mourisl/centrifuger) and run `make` in it:

```bash
git clone https://github.com/mourisl/centrifuger.git
cd centrifuger
make
```

The executables are written into the directory you just built in. There is no `make install` step.

### Putting Centrifuger on your PATH

After building, you can run the programs by their full path — `/path/to/centrifuger/centrifuger` — but it is more convenient to make them available everywhere. Either add the build directory to `PATH`:

```bash
export PATH="/path/to/centrifuger:$PATH"
```

(add that line to your `~/.bashrc` or `~/.zshrc` to make it permanent), or create a symbolic link from a directory that is already on `PATH`:

```bash
ln -s /path/to/centrifuger/centrifuger ~/bin/centrifuger
```

:::caution Link every program you need
`centrifuger` alone is not enough for a full workflow. Link or expose `centrifuger-build`, `centrifuger-quant` and `centrifuger-download` as well.
:::

## Verify the installation

Run a program without arguments to print its usage message:

```bash
centrifuger
centrifuger-build
centrifuger-quant
```

For an end-to-end check, use the `example/` directory included in the source distribution — it contains a tiny reference, a taxonomy and a pair of FASTQ files, and the whole test takes seconds. The [quick start](/getting-started/quick-start/#run-the-bundled-example) walks through it.

## Keeping up to date

Conda installations update in the usual way:

```bash
conda update -c conda-forge -c bioconda centrifuger
```

Source installations update by pulling and rebuilding:

```bash
cd /path/to/centrifuger
git pull
make clean && make
```

:::note Indexes and versions
Indexes are not rebuilt by a software update, and Centrifuger indexes are portable across machines. If a release ever changes the index format, the release notes on GitHub will say so; otherwise an existing `*.cfr` index keeps working after an upgrade.
:::

## Next steps

- [Quick start](/getting-started/quick-start/) — your first classification.
- [Pre-built indexes](/guides/prebuilt-indexes/) — download a database instead of building one.
