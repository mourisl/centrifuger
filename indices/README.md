### make
Combined pipeline for centrifuger-download and centriduger-build to download and create indices for sequences from NCBI RefSeq or GenBank.

### SILVA
run `perl download-silva.pl -h` for help message and default parameter.

run `perl download-silva.pl` to obtain the nodes.dmp, names.dmp, silva_seqid_to_taxid.map, and silva_seq.fa.gz files for the input of centrifuger-build.

### GTDB
run `perl download-gtdb.pl` to obtain the files to create GTDB index. The files are gtdb_nodes.dmp, gtdb_names.dmp, gtdb_file.list, gtdb_seqid2taxid.map and the folder of GTDB's representative genome files.

If you have NCBI's names.dmp file, e.g. from "centrifuger-download taxonomy", you can add the "--names names.dmp", which will try to match the taxonomy ID for intermediate taxonomy levels.
