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
A pre-built index for GTDB r226 representative genomes is available at: https://www.dropbox.com/scl/fo/xjp5r81jxkzxest9ijxul/ADfYFKoxIyl0hrICeEI63QM?rlkey=5lij0ocrbre165pa52mavux5z&st=4ol28yv2&dl=0 . 

### Create seqid2taxid.map from NCBI accession2taxid mapping file
The seqID, e.g. accession ID, to taxonomy ID mapping file is necessary for Centrifuge and optional for Centrifuger. In case you have downloaded the genomes and also have access to the (nucl_gb.accession2taxid.gz)[https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz] file, here is the procedure to create the seqid_taxid.map file.

```
find genomes -type f -name "*.fna.gz" | xargs -I{} zcat {} | grep "^>" | cut -f1 -d' ' | tr -d ">" > seqid.list # Get the seqID list

perl SearchAccessionIdToTaxId.pl seqid.list <(zcat nucl_gb.accession2taxid.gz) > seqid_to_taxid.map # seqIDs not in the accession file will be assigned to taxonomy ID 1 in this script.
```

In the case you want to filter the genomes that is not in the accession2taxid mapping file, you can use the following steps:
```
find genomes -type f -name "*.fna.gz" | xargs -I{} zgrep -H "^>" {} | cut -f1 -d' ' | tr -d ">" | tr ":" "\t" > seqid_wfname.list # Create a two-column file where the first column is the file name and the second column is the seqID.
perl SearchAccessionIdToTaxId.pl <(cut -f2 seqid_wfname.list) <(zcat nucl_gb.accession2taxid.gz) > seqid_to_taxid.map
paste seqid_wfname.list seqid_to_taxid.map | awk '{if ($4!=1) print $1,$4}' | uniq > genome_file_wtaxid.list 
```
The genomes_file_wtaxid.list file can be used for the "-l" option in centrifuger-build and bypass the requirement of the "--conversion-table". Using awk on "$2,$4" will create the filtered seqid_to_taxid.map file.
