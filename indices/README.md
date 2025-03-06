### make
Combined pipeline for centrifuger-download and centriduger-build to download and create indices for sequences from NCBI RefSeq or GenBank.

### SILVA
run `perl silva-download.pl -h` for help message and default parameter.

run `perl silva-download.pl` to obtain the nodes.dmp, names.dmp, silva_seqid_to_taxid.map, and silva_seq.fa.gz files for the input of centrifuger-build.

### GTDB
run `perl gtdb-download.pl` to obtain the files to create GTDB index. The files are gtdb_nodes.dmp, gtdb_names.dmp, gtdb_file.list, gtdb_fname_to_taxid.map and the folder of GTDB's representative genome files.

If you have NCBI's names.dmp file, e.g. from "centrifuger-download taxonomy", you can add the "--names names.dmp", which will try to match the taxonomy ID for intermediate taxonomy levels. If you need to the sequence ID to taxonomy ID mapping file, e.g. for Centrifuge's index, you can add the "--generateSeqId2TaxId" option.

The repetitiveness of GTDB representative genomes is very low, so the run-block compression does not provide much space saving. Therefore, it is recommended to add the "--rbbwt-b 1" option in centrifuger-build to disable the compression, which can accelerate the classification procedure. Due to the large volume of the database, you can increase the SA sampling sparsity (--offrate option) to reduce the index size. For example, you can run centrifuger-build as

```
../centrifuger-build -l gtdb_fname_to_taxid.map --taxonomy-tree gtdb_nodes.dmp --name-table gtdb_names.dmp -t 16 -o cfr_gtdb --build-mem 500G --offrate 5 --rbbwt-b 1
```
A pre-built index for GTDB r220 representative genomes is available at: https://www.dropbox.com/scl/fo/wtsermtdd62n1ttryuhwv/ADWTYiJA6Dh5gdbHmwC2nfo?rlkey=fgdpthjukcj1uhsjrcv35j927&dl=0 .

### Create seqid2taxid.map from NCBI accession2taxid mapping file
The seqID, e.g. accession ID, to taxonomy ID mapping file is necessary for Centrifuge and optional for Centrifuger. In case you have downloaded the genomes and also have access to the (nucl_gb.accession2taxid.gz)[https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz] file, here is the procedure to create the seqid_taxid.map file.

```

find genomes -type f -name "*.fna.gz" | xargs -I{} zcat {} | grep "^>" | cut -f1 -d' ' | tr -d ">" > seqid.list # Get the seqID list

perl SearchAccessionIdToTaxId.pl seqid.list <(zcat nucl_gb.accession2taxid.gz) > seqid_to_taxid.map # seqIDs not in the accession file will be assigned to taxonomy ID 1 in this script.
```
