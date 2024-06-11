Centrifuger
=======

Described in: 

  Song, L., Langmead B.. Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification. Genome Biol. 2024 Apr 25;25(1):106. doi: 10.1186/s13059-024-03244-4.

  Copyright (C) 2023-, Li Song


### What is Centrifuger?

Centrifuger is an efficient taxonomic classification method that compares sequencing reads against a microbial genome database. It implemented a novel lossless compression method, run-block comprssed BWT, and other strategies to efficiently reduce the size of the microbial genome database like RefSeq. For example, Centrifuger can classify reads against the 2023 RefSeq prokaryotic genomes containing about 140G nucleotides using 43 GB memory. Despite running on a compressed data structure, Centrifuger is also highly efficient and can process a typical sequencing sample within an hour. 

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/centrifuger), e.g. with `git clone https://github.com/mourisl/centrifuger.git`
2. Run `make` in the repo directory

You will find the executable files in the downloaded directory. If you want to run Centrifuger without specifying the directory, you can either add the directory of Centrifuger to the environment variable PATH or create a soft link ("ln -s") of the file "centrifuger" to a directory in PATH.

Centrifuger depends on [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads). 

Centrifuger is also available from [Bioconda](https://anaconda.org/bioconda/centrifuger). You can install Centrifuger with `conda install -c conda-forge -c bioconda centrifuger`.

### Usage

#### Build index

    Usage: ./centrifuger-build [OPTIONS]
      Required:
        -r FILE: reference sequence file (can use multiple -r to specify more than one input file)
            or
        -l FILE: list of reference sequence file stored in <file>, one sequence file per row. Can include the taxonomy ID mapping information in the second column.
        --taxonomy-tree FILE: taxonomy tree, i.e., nodes.dmp file
        --name-table FILE: name table, i.e., names.dmp file
      Optional:
        --conversion-table FILE: seqID to taxID conversion file
          When not set, expect -l option and the -l file should have two columns as "file taxID"
        -o STRING: output prefix [centrifuger]
        -t INT: number of threads [1]
        --build-mem STR: automatic infer bmax and dcv to match memory constraints, can use T,G,M,K to specify the memory size [not used]
        --bmax INT: block size for blockwise suffix array sorting [16777216]
        --dcv INT: difference cover period [4096]
        --offrate INT: SA/offset is sampled every (2^<int>) BWT chars [4]
        --subset-tax INT: only consider the subset of input genomes under taxonomy node <int> [0]

An example of pre-built index containing human, bacteria, archea, and virus genomes from RefSeq plus SARS-CoV-2 variants from GenBank is available at [Zenodo](https://zenodo.org/records/10023239). The default --bmax and --dcv option may be inefficient for building indexes for larger genome databases, please use --build-mem option to specify the rough estimation of the available memory.

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
        --barcode STR: path to the barcode file
        --UMI STR: path to the UMI file
        --read-format STR: format for read, barcode and UMI files, e.g. r1:0:-1,r2:0:-1,bc:0:15,um:16:-1 for paired-end files with barcode and UMI
        --min-hitlen INT: minimum length of partial hits [auto]
        --hitk-factor INT: resolve at most <int>*k entries for each hit [40; use 0 for no restriction]
        --merge-readpair: merge overlapped paired-end reads and trim adapters 

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
The sixth column is the number of base pairs of the read that match the genomic sequence found by Centrifuger (e.g., 80) 
The seventh column is the length of a read or the combined length of mate pairs (e.g., 80). 
The eighth column is the number of classifications for this read, indicating how many assignments were made in the output (e.g.,1).
```

### Practical notes
* #### Create index for genomes from NCBI.

You can use "centrifuger-download" to download reference sequences from NCBI. The following two commands download the NCBI taxonomy to taxonomy/ in the current directory, and all complete archaeal, bacterial and viral genomes to library/.

	./centrifuger-download -o taxonomy taxonomy
	./centrifuger-download -o library -d "archaea,bacteria,viral" refseq > seqid2taxid.map

To add human (taxonomy ID 9606) or mouse (taxonomy ID 10090) genome to the downloaded files, you can use the following command

	# human: T2T-CHM13
	./centrifuger-download -o library -d "vertebrate_mammalian" -t 9606 refseq >> seqid.map
	# human: hg38 reference genome
	./centrifuger-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' refseq	
	# mouse
	./centrifuger-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 10090 -c 'reference genome' refseq >> seqid.map

To build the index, first put the downloaded files in a list (this part is different from Centrifuge, where the files need to be concatendated) and then run centrifuger-build:
	
	ls library/*/*.fna.gz > file.list # use *_dustmasked.fna.gz as the file list if using dustmasker in centrifuger-download 

	## build centrifuger index with 4 threads on a server with 300GB memory
	./centrifuger-build -t 4 --conversion-table seqid2taxid.map \
		--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		-l file.list -o refseq_abv --build-mem 240G
	
After building the index, all but the refseq_abv.[1234].cfr index files may be removed.

* #### Build custom database index 
The index building procedure is similar to [Centrifuge's](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building), but with names changing to centrifuger. For example, centrifuge-download is centrifuger-download. 

* #### 10x Genomics data and barcode-based single-cell data

If your input has barcode information, you can use "--barcode" to specify the barcode file and use "--read-format" to tell Centrifuger how to extract barcode information. The "--read-format" option can also specify the extraction for read1, read2 and UMI. The value for this argument is a comma-separated string, each field in the string is also a semi-comma-splitted string

	[r1|r2|bc|um]:start:end:strand

The start and end are inclusive and -1 means the end of the read. You may use multiple fields to specify non-consecutive segments, e.g. bc:0:15,bc:32:-1. The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction. The strand symbol can be omitted if it is '+' and is ignored on r1 and r2. For example, when the barcode is in the first 16bp of read1, one can use the option `-1 read1.fq.gz -2 read2.fq.gz --barcode read1.fq.gz --read-format bc:0:15,r1:16:-1`.

Centrifuger supports using wildcard in the -1 -2/-u option, so a typical way to run 10x Genomics single-end data is by:

	./centrifuger -x cfr_idx -u "path_to_10x_fastqs/*_R2_*.fastq.gz" --barcode "path_to_10x_fastqs/*_R1_*.fastq.gz" --UMI "path_10x_fastqs/*_R1_*.fastq.gz" --read-format bc:0:15,um:16:-1 --barcode-whitelist cellranger_folder/cellranger-cs/VERSION/lib/python/cellranger/barcodes/3M-february-2018.txt.gz [other options]

The exact options depend on your 10x Genomics kit. The quotes around the paths with wildcard  are necessary.

Moreover, Centrifuger can translate input cell barcodes to another set of barcodes. You can specify the translation file through the option --barcodeTranslate. The translation file is a two-column tsv/csv file with the translated barcode on the first column and the original barcode on the second column. This option also supports combinatorial barcoding, such as SHARE-seq. Centrifuger can translate each barcode segment provided in the second column to the ID in the first column and add "-" to concatenate the IDs in the output.

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
