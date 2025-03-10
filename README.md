Centrifuger
=======

Described in: 

  Song, L., Langmead B.. Centrifuger: lossless compression of microbial genomes for efficient and accurate metagenomic sequence classification. Genome Biol. 2024 Apr 25;25(1):106. doi: 10.1186/s13059-024-03244-4. **Best Paper Award at RECOMB2024**

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

An example of pre-built index containing human, bacteria, archea, and virus genomes from RefSeq plus SARS-CoV-2 variants from GenBank is available at [Zenodo](https://zenodo.org/records/10023239). Other pre-built indices are [GTDB_r220](https://gtdb.ecogenomic.org/) available at [Dropbox](https://www.dropbox.com/scl/fo/wtsermtdd62n1ttryuhwv/ADWTYiJA6Dh5gdbHmwC2nfo?rlkey=fgdpthjukcj1uhsjrcv35j927&dl=0), and NCBI nt database from 2018 available at [Dropbox](https://www.dropbox.com/scl/fo/dbjzuerib7nfluj2u3yw7/AER-MP9RaqL0g59YwTXe4RU?rlkey=7zy78nvwdw19fuaetnn3u7qhi&st=6un1f4ux&dl=0). The default --bmax and --dcv option may be inefficient for building indexes for larger genome databases, please use --build-mem option to specify the rough estimation of the available memory.

#### Classification

    Usage: ./centrifuger [OPTIONS] > classification.tsv
      Required:
        -x FILE: index prefix
        -1 FILE -2 FILE: paired-end read files
          or
        -u FILE: single-end read file
          or
        -i FILE: interleaved paried-end read file
      Optional:
        -t INT: number of threads [1]
        -k INT: report upto <int> distinct, primary assignments for each read pair [1]
        --un STR: output unclassified reads to files with the prefix of <str>, e.g. <str>_1/2.fq.gz
        --cl STR: output classified reads to files with the prefix of <str>
        --barcode STR: path to the barcode file
        --UMI STR: path to the UMI file
        --read-format STR: format for read, barcode and UMI files, e.g. r1:0:-1,r2:0:-1,bc:0:15,um:16:-1 for paired-end files with barcode and UMI
        --min-hitlen INT: minimum length of partial hits [auto]
        --hitk-factor INT: resolve at most <int>*k entries for each hit [40; use 0 for no restriction]
        --merge-readpair: merge overlapped paired-end reads and trim adapters 

#### Quantification (taxonomic profiling)

    Usage: ./centrifuger-quant [OPTIONS] > report.tsv
      Required:
        -x FILE: index prefix
        -c FILE: classification file
      optional:
        --min-score INT: only consider reads with score at least <int> 
        --min-length INT: only consider reads with classified length at least <int>
        --output-format INT: output format. (0:centrifuge,default, 1:Metaphlan, 2:CAMI)        

The quantification results will be affected by the "-k" option from the classification program "centrifuger". Increasing "-k" will provide ambiguous but more specific classification result, potentially can improve the quantification result.  

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

The "centrifuger-quant" estimate the abundance for each taxonomy ID, and the quantification output has 7 columns.

```
name	taxID	taxRank	genomeSize	numReads	numUniqueReads	abundance
Legionella_pneumophila_subsp._pneumophila_str._Philadelphia_1	272624	strain	3397753	50	48	0.392641

The first column is the name of a genome, or the name corresponding to a taxonomic ID (the second column) at a rank higher than the strain (e.g., Legionella_pneumophila_str._Pari).
The second column is the taxonomic ID (e.g., 297246).
The third column is the taxonomic rank (e.g., strain).
The fourth column is the length of the genome sequence (e.g., 3503503).
The fifth column is the number of reads classified to some genomic sequences (multi-classified reads are evenly distributed) under this taxonomy node (e.g., 50).
The sixth column is the number of reads uniquely classified to a genomic sequence under this taxonomy node (e.g., 48).
The seventh column is the proportion of this genome normalized by its genomic length (e.g., 0.392641).
```

### Practical notes
* #### Create index for genomes from NCBI.

You can use "centrifuger-download" to download reference sequences from NCBI. The following two commands download the NCBI taxonomy to taxonomy/ in the current directory, and all complete archaeal, bacterial and viral genomes to library/.

	./centrifuger-download -o taxonomy taxonomy
	./centrifuger-download -o library -d "archaea,bacteria,viral" refseq > seqid2taxid.map

To add human (taxonomy ID 9606) or mouse (taxonomy ID 10090) genome to the downloaded files, you can use the following command

	# human: T2T-CHM13
	./centrifuger-download -o library -d "vertebrate_mammalian" -t 9606 refseq >> seqid2taxid.map
	# human: hg38 reference genome
	./centrifuger-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' refseq >> seqid2taxid.map
	# mouse
	./centrifuger-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 10090 -c 'reference genome' refseq >> seqid2taxid.map

To build the index, first put the downloaded files in a list (this part is different from Centrifuge, where the files need to be concatendated) and then run centrifuger-build:
	
	find library -type f -name "*.fna.gz" > file.list # use *_dustmasked.fna.gz as the file list if using dustmasker in centrifuger-download 

	## build centrifuger index with 4 threads on a server with 300GB memory
	./centrifuger-build -t 4 --conversion-table seqid2taxid.map \
		--taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
		-l file.list -o refseq_abv --build-mem 240G
	
After building the index, all but the refseq_abv.[1234].cfr index files may be removed.

* #### Build custom database index 
The index building procedure is similar to [Centrifuge's](http://www.ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building), but with names changing to centrifuger. For example, centrifuge-download is centrifuger-download. 

The folder "indices" contains information for creating index from other sources, like SILVA and GTDB.

* #### 10x Genomics data and barcode-based single-cell data

If your input has barcode information, you can use "--barcode" to specify the barcode file and use "--read-format" to tell Centrifuger how to extract barcode information. The "--read-format" option can also specify the extraction for read1, read2 and UMI. The value for this argument is a comma-separated string, each field in the string is also a semi-comma-splitted string

	[r1|r2|bc|um]:start:end:strand

The start and end are inclusive and -1 means the end of the read. You may use multiple fields to specify non-consecutive segments, e.g. bc:0:15,bc:32:-1. The strand is presented by '+' and '-' symbol, if '-' the barcode will be reverse-complemented after extraction. The strand symbol can be omitted if it is '+' and is ignored on r1 and r2. For example, when the barcode is in the first 16bp of read1, one can use the option `-1 read1.fq.gz -2 read2.fq.gz --barcode read1.fq.gz --read-format bc:0:15,r1:16:-1`. If "--barcode" or "--UMI" option is ignored in the command, Centrifuger will extract the pattern from read1. Note that "--barcode-whilelist" option requires the use of "--barcode" option to specify the barcode file.

Centrifuger supports using wildcard in the -1 -2/-u option, so a typical way to run 10x Genomics single-end data is by:

	./centrifuger -x cfr_idx -u "path_to_10x_fastqs/*_R2_*.fastq.gz" \
		--barcode "path_to_10x_fastqs/*_R1_*.fastq.gz" --UMI "path_10x_fastqs/*_R1_*.fastq.gz" 
		--read-format bc:0:15,um:16:-1 \
		--barcode-whitelist cellranger_folder/cellranger-cs/VERSION/lib/python/cellranger/barcodes/3M-february-2018.txt.gz [other options]

The exact options depend on your 10x Genomics kit. The quotes around the paths with wildcard  are necessary.

Moreover, Centrifuger can translate input cell barcodes to another set of barcodes. You can specify the translation file through the option --barcodeTranslate. The translation file is a two-column tsv/csv file with the translated barcode on the first column and the original barcode on the second column. This option also supports combinatorial barcoding, such as SHARE-seq. Centrifuger can translate each barcode segment provided in the second column to the ID in the first column and add "-" to concatenate the IDs in the output.

The bc and um option can parse the barcode and UMI from the fastq header comment field. The format is [bc|um]:hd:field:start:end:strand. "hd" is a keyword so the search will be in the header comment. "field" can be a number (0-based), which is specifies which field in the comment (read id is excluded) contains the barcode/UMI. "field" can also be a string, and it search for the pattern starting with the "field" and extract the barcode/UMI from there. For example, if the header looks like "@r1 CR:Z:NNNN CB:Z:ACGT UR:Z:NNNN", then "bc:hd:1:5:-1" or "bc:hd:CB:5:-1" will extract the barcode "ACGT" from the header. 

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

### centrifuger classification command examples
 
```
  ./centrifuger -x cfr_idx -u read_single_end.fq > output.tsv # single-end data
  ./centrifuger -x cfr_idx -1 read_1.fq -2 read_2.fq > output.tsv # paired-end data
  ./centrifuger -x cfr_idx -u "*_R2_*.fastq.gz" --barcode "*_R1_*.fastq.gz" --UMI "*_R1_*.fastq.gz" --read-format bc:0:15,um:16:-1 --barcode-whitelist 3M-february-2018.txt.gz > output.tsv # 10x Genomics single-end fastq file
  samtools fastq -T CB -f 4 alignment.bam | ./centrifuger -x cfr_idx --read-format bc:hd:CB:5:-1 -u - > output.tsv # unmapped reads from single-end alignment BAM file generated by 10x Genomics Cellranger.
  samtools view -f 0xC alignment.bam | samtools sort -n - | samtools fastq -T CB,UB - | ./centrifuger -x cfr_idx --read-format bc:hd:CB:5:-1,um:hd:UB:5:-1 --un unclassified -i - > output.tsv # unmapped read fragment from paired-end alignment BAM file generated by 10x Genomics Cellranger, including barcode and UMI in the output. Also output unclassified reads to files with prefix unclassified.
```

### Support

Create a [GitHub issue](https://github.com/mourisl/centrifuger/issues). 
