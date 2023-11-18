Centrifuger
=======

Described in: 

  Song, L., Langmead B.. Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification. Biorxiv (2023). https://doi.org/10.1101/2023.11.15.567129

  Copyright (C) 2023-, Li Song


### What is Centrifuger?

Centrifuger is an efficient taxonomic classification method that compares sequencing reads against a microbial genome database. It implemented a novel lossless compression method, run-block comprssed BWT, and other strategies to efficiently reduce the size of the microbial genome database like RefSeq. For example, Centrifuger can classify reads against the 2023 RefSeq prokaryotic genomes containing about 140G nucleotides using 43 GB memory. Despite running on a compressed data structure, Centrifuger is also highly efficient and can process a typical sequencing sample within an hour. 

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/centrifuger), e.g. with `git clone --recurse-submodules https://github.com/mourisl/centrifuger.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run Centrifuger without specifying the directory, you can either add the directory of Centrifuger to the environment variable PATH or create a soft link ("ln -s") of the file "centrifuger" to a directory in PATH.

Centrifuger depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads). 

### Usage

#### Build index

    Usage: ./centrifuger-build [OPTIONS]
      Required:
        -r FILE: reference sequence file (can use multiple -r to specify more than one input file)
            or
        -l FILE: list of reference sequence file stored in <file>, one sequence file per row
        --taxonomy-tree FILE: taxonomy tree, i.e., nodes.dmp file
        --name-table FILE: name table, i.e., names.dmp file
        --conversion-table FILE: seqID to taxID conversion file
      Optional:
        -o STRING: output prefix [centrifuger]
        -t INT: number of threads [1]
        --bmax INT: block size for blockwise suffix array sorting [16777216]
        --offrate INT: SA/offset is sampled every (2^<int>) BWT chars [4]
        --dcv INT: difference cover period [4096]
        --build-mem STR: automatic infer bmax and dcv to match memory constraints, can use P,G,M,K to specify the memory size [not used]
        --subset-tax INT: only consider the subset of input genomes under taxonomy node <int> [0]

An example of pre-built index containing human, bacteria, archea, and virus genomes from RefSeq plus SARS-CoV-2 variants from GenBank is available at [Zenodo](https://zenodo.org/records/10023239).

#### Classification

    Usage: ./centrifuger [OPTIONS]
      Required:
        -x FILE: index prefix
        -1 FILE -2 FILE: paired-end read
          or
        -u FILE: single-end read
      Optional:
        -o STRING: output prefix [centrifuger]
        -t INT: number of threads [1]
        -k INT: report upto <int> distinct, primary assignments for each read pair [1]
        --min-hitlen INT: minimum length of partial hits [auto]
        --hitk-factor INT: resolve at most <int>*k entries for each hit [40; use 0 for no restriction]

### Input/Output

The primary input to Centrifuger is the index of the genome database (-x), and gzipped or uncompressed read fastq files (-1/-2 for paired; -u for single-end).

The output is to stdout, with the TSV format as following:
```
readID    seqID   taxID score      2ndBestScore    hitLength    queryLength numMatches
1_1       MT019531.1     2697049   4225       0               80   80      1

The first column is the read ID from a raw sequencing read (e.g., 1_1 in the example).
The second column is the sequence ID of the genomic sequence, where the read is classified (e.g., MT019531.1).
The third column is the taxonomic ID of the genomic sequence in the second column (e.g., 2697049).
The fourth column is the score for the classification, which is the weighted sum of hits (e.g., 4225)
The fifth column is the score for the next best classification (e.g., 0).
The sixth column is a pair of two numbers: (1) an approximate number of base pairs of the read that match the genomic sequence and (2) the length of a read or the combined length of mate pairs (e.g., 80 / 80).
The seventh column is a pair of two numbers: (1) an approximate number of base pairs of the read that match the genomic sequence and (2) the length of a read or the combined length of mate pairs (e.g., 80 / 80). 
The eighth column is the number of classifications for this read, indicating how many assignments were made (e.g.,1).
```

### Practical notes
#### Build custom database index 
The index building procedure is similar to [Centrifuge's](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building), but with names changing to centrifuger. For example, centrifuge-download is centrifuger-download. 

### Example

The directory "./example" in this distribution contains files for building Centrifuger index and classification. Suppose you are in the example folder, and Centrifuger has been compiled with "make" command.

* Build index:
```
../centrifuger-build -r ref.fa --taxonomy-tree nodes.dmp --name-table names.dmp --conversion-table ref_seqid.map -o cfr_ref_idx
```  
After running the above command, you shall see the index file "cfr_ref_idx.*.cfr" in the example folder.

* Classification
```
../centrifuger -1 example_1.fq -2 example_2.fq -x cfr_ref_idx > output.tsv
```
The output.tsv should be similar to the example_class.out file in the folder. 

### Support

Create a [GitHub issue](https://github.com/mourisl/centrifuger/issues). 
