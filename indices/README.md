### make
Combined pipeline for centrifuger-download and centriduger-build to download and create indices for sequences from NCBI RefSeq or GenBank.

### SILVA
run `perl download-silva.pl -h` for help message and default parameter.

run `perl download-silva.pl` to obtain the nodes.dmp, names.dmp, silva_seqid_to_taxid.map, and silva_seq.fa.gz files for the input of centrifuger-build.

### GTDB
run `perl download-gtdb.pl` to obtain the files to create GTDB index. The files are gtdb_nodes.dmp, gtdb_names.dmp, gtdb_file.list, gtdb_seqid2taxid.map and the folder of GTDB's representative genome files.

If you have NCBI's names.dmp file, e.g. from "centrifuger-download taxonomy", you can add the "--names names.dmp", which will try to match the taxonomy ID for intermediate taxonomy levels.

The repetitiveness of GTDB representative genomes is very low, so the run-block compression does not provide much space saving. Therefore, it is recommended to add the "--rbbwt-b 1" option in centrifuger-build to disable the compression, which can accelerate the classification procedure. Due to the large volume of the database, you can increase the SA sampling sparsity (--offrate option) to reduce the index size. For example, you can run centrifuger-build as

```
../centrifuger-build -l gtdb_file.list --taxonomy-tree ./gtdb_nodes.dmp --name-table ./gtdb_names.dmp --conversion-table gtdb_seqid_to_taxid.map -t 16 -o cfr_gtdb --build-mem 500G --offrate 5 --rbbwt-b 1
```
