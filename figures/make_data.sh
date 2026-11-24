#!/bin/bash
# Input tables for the README figures, from the example data in this repository.
# usage: figures/make_data.sh UPSTREAM_QUANT
#   UPSTREAM_QUANT: centrifuger-quant built from mourisl/centrifuger v1.2.0, for the "before" counts.
# Writes figures/shared_reads/counts.tsv and figures/protein_weights/weights.tsv.
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
X=$R/example
UP=$1
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
make -C "$R" -s centrifuger centrifuger-build centrifuger-quant
Q=$R/centrifuger-quant

# Figure 1: one read that fits three sibling strains (genome lengths 1, 1 and 2 Mb), as in tests/run_quant_tests.sh.
cp "$X/nodes.dmp" "$T/nodes3.dmp"; cp "$X/names.dmp" "$T/names3.dmp"
printf '900001\t|\t446\t|\tstrain\t|\n' >> "$T/nodes3.dmp"
printf '900001\t|\tStrain_C\t|\t\t|\tscientific name\t|\n' >> "$T/names3.dmp"
printf '272624 1000000\n297246 1000000\n900001 2000000\n' > "$T/sizes3.tsv"
{
  printf 'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
  for t in 272624 297246 900001; do printf 'r1\ts\t%d\t100\t100\t100\t100\t3\n' "$t"; done
} > "$T/shared.tsv"
TX=(--taxonomy-tree "$T/nodes3.dmp" --name-table "$T/names3.dmp" --size-table "$T/sizes3.tsv" -c "$T/shared.tsv")
"$Q" "${TX[@]}" > "$T/f0.tsv" 2>/dev/null
"$UP" "${TX[@]}" --output-format 3 > "$T/up.f3.tsv" 2>/dev/null
"$Q" "${TX[@]}" --output-format 3 > "$T/new.f3.tsv" 2>/dev/null
"$Q" "${TX[@]}" --output-format 3 --em-reads > "$T/em.f3.tsv" 2>/dev/null
mkdir -p "$R/figures/shared_reads"
awk -F'\t' -v OFS='\t' -v F0="$T/f0.tsv" -v UPF="$T/up.f3.tsv" -v NEWF="$T/new.f3.tsv" -v EMF="$T/em.f3.tsv" '
  function kr(file, key,   line, f) { while ((getline line < file) > 0) { split(line, f, "\t"); v[key, f[5]] = f[3]; s[key] += f[3] } close(file) }
  BEGIN {
    while ((getline line < F0) > 0) { split(line, f, "\t"); if (f[2] ~ /^[0-9]+$/) { v["numReads", f[2]] = f[5]; v["emNumReads", f[2]] = f[8] } }
    kr(UPF, "kreport_v1.2.0"); kr(NEWF, "kreport"); kr(EMF, "kreport_em_reads")
    s["numReads"] = 1; s["emNumReads"] = 1   # root clade reads in format 0
    print "output", "Strain_A_1Mb", "Strain_B_1Mb", "Strain_C_2Mb", "taxon_reads_sum"
    n = split("numReads emNumReads kreport_v1.2.0 kreport kreport_em_reads", k, " ")
    for (i = 1; i <= n; i++) print k[i], v[k[i], 272624] + 0, v[k[i], 297246] + 0, v[k[i], 900001] + 0, s[k[i]] + 0
  }' > "$R/figures/shared_reads/counts.tsv"

# Figure 2: the example classified against a protein-marked copy of the example index (hit lengths in amino acids),
# read as nucleotide (v1.2.0 behaviour, and the fork without --protein) and as protein.
"$R/centrifuger-build" -r "$X/ref.fa" --taxonomy-tree "$X/nodes.dmp" --name-table "$X/names.dmp" \
  --conversion-table "$X/ref_seqid.map" -o "$T/ex" > /dev/null 2>&1
"$R/centrifuger" -x "$T/ex" -1 "$X/example_1.fq" -2 "$X/example_2.fq" -k 5 > "$T/class.tsv" 2> /dev/null
for i in 1 2 3; do cp "$T/ex.$i.cfr" "$T/aa.$i.cfr"; done
awk 'BEGIN{FS=OFS="\t"} $1=="sequence_type"{$2="amino_acid"} 1' "$T/ex.4.cfr" > "$T/aa.4.cfr"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {h=int($6/3); if(h>66)h=66; $6=h; print}' "$T/class.tsv" > "$T/class_aa.tsv"
"$R/centrifuger-inspect" -x "$T/ex" --taxonomy-tree > "$T/nodes.dmp"
"$R/centrifuger-inspect" -x "$T/ex" --name-table > "$T/names.dmp"
"$R/centrifuger-inspect" -x "$T/ex" --size-table > "$T/size.tsv"
TA=(--taxonomy-tree "$T/nodes.dmp" --name-table "$T/names.dmp" --size-table "$T/size.tsv" -c "$T/class_aa.tsv")
"$Q" -x "$T/ex" -c "$T/class.tsv" > "$T/nt.tsv" 2>/dev/null
"$Q" "${TA[@]}" > "$T/aa_as_nt.tsv" 2>/dev/null
"$UP" "${TA[@]}" > "$T/aa_v1.2.0.tsv" 2>/dev/null
"$Q" "${TA[@]}" --protein > "$T/aa_protein.tsv" 2>/dev/null
"$Q" -x "$T/aa" -c "$T/class_aa.tsv" > "$T/aa_index.tsv" 2>/dev/null
mkdir -p "$R/figures/protein_weights"
{
  printf 'run\ttaxID\tname\tnumReads\temNumReads\temNumReadsWeighted\tabundance\n'
  for run in nt aa_v1.2.0 aa_as_nt aa_protein aa_index; do
    awk -F'\t' -v OFS='\t' -v run=$run 'NR>1 && $2 ~ /^[0-9]+$/ && ($2==272624 || $2==297246 || $2==1) {print run, $2, $1, $5, $8, $9, $7}' "$T/$run.tsv"
  done
} > "$R/figures/protein_weights/weights.tsv"
cat "$R/figures/shared_reads/counts.tsv" "$R/figures/protein_weights/weights.tsv"
