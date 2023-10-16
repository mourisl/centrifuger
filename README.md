Centrifuger
=======

Described in: 

  Song, L., Langmead B.. Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification.

  Copyright (C) 2023-, Li Song


### What is Centrifuger?

Centrifuger is a 

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/centrifuger), e.g. with `git clone https://github.com/mourisl/centrifuger.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run Centrifuger without specifying the directory, you can either add the directory of Centrifuger to the environment variable PATH or create a soft link ("ln -s") of the file "centrifuger-class" to a directory in PATH.

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

#### Classification

  Usage: ./centrifuger-class [OPTIONS]
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
    

### Practical notes
#### Download your own sequences
The download procedure is similar to [Centrifuge's](), but with names changing to centrifuger. For example, centrifuge-download is centrifuger-download. 

### Example



### Support

Create a [GitHub issue](https://github.com/mourisl/centrifuger/issues). 
