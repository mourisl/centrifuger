## taxonomy.py
A python script to obtain various information from the nodes.dmp file downloaded from NBCI. Run "python3 taxonomy.py -h" for more information

## requant-centrifuge.pl
The perl script that can use centrifuger-quant to recalculate the abundance from Centrifuge classification result. 
It takes four parameters, in the order of: path to Centrifuge, path to Centrifuger, Centrifuge's index prefix, Centrifuge's classification output.
An example is:
```
  perl requant-centrifuge.pl ./centrifuge/ ./centrifuger/ cf_idx/p_compressed+h+v cf_classification.out > new_report.tsv
```
The idea of the script is to run "centrifuge-inspect" to get the taxonomy tree, name, and genome size from Centrifuge's index, and then run "centrifuger-quant" to get the report.
