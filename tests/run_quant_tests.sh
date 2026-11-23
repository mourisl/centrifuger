#!/bin/bash
# Tests for centrifuger-quant: Quantifier state under AddressSanitizer, kraken-report read counts that add up (with
# and without --em-reads), and read weights for a protein index.
# usage: tests/run_quant_tests.sh            (needs make, g++ with -fsanitize=address, zlib)
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
X=$R/example
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fail=0
ok() { printf '%-66s ok\n' "$1" ; }
bad() { printf '%-66s FAIL\n' "$1" ; fail=1 ; }

make -C "$R" -s centrifuger centrifuger-build centrifuger-quant
Q=$R/centrifuger-quant

# Example index and classification (as in the README), keeping up to 5 hits per read.
"$R/centrifuger-build" -r "$X/ref.fa" --taxonomy-tree "$X/nodes.dmp" --name-table "$X/names.dmp" \
  --conversion-table "$X/ref_seqid.map" -o "$T/ex" > "$T/build.log" 2>&1
"$R/centrifuger" -x "$T/ex" -1 "$X/example_1.fq" -2 "$X/example_2.fq" -k 5 > "$T/class.tsv" 2> "$T/class.log"

# The same index marked as protein, and the classification with hit lengths in amino acids
# (a 2x100 nt pair translates to at most 2x33 amino acids).
for i in 1 2 3; do cp "$T/ex.$i.cfr" "$T/aa.$i.cfr"; done
sed 's/^sequence_type\tnucleotide/sequence_type\tamino_acid/' "$T/ex.4.cfr" > "$T/aa.4.cfr"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {h=int($6/3); if(h>66)h=66; $6=h; print}' "$T/class.tsv" > "$T/class_aa.tsv"

echo "== Quantifier state (AddressSanitizer)"
g++ -std=c++11 -g -O1 -fsanitize=address -fno-omit-frame-pointer -Wno-unused-result -I"$R" "$R/tests/quant_state_test.cpp" \
  -o "$T/quant_state_test" -lz -lpthread
"$T/quant_state_test" "$T/ex" "$T/aa" "$T/class.tsv" "$X/nodes.dmp" "$X/names.dmp" 2> "$T/asan.log" || { fail=1; cat "$T/asan.log"; }

# Every row's clade reads are its taxon reads plus its children's clade reads, and the taxon reads add up to the
# root's clade reads (Pavian divides by their sum). Depth comes from the name's indentation, two spaces a level.
kreport_adds_up() {
  awk -F'\t' '
    { match($6, /^ */); dep = RLENGTH / 2; clade[NR] = $2; taxon[NR] = $3; childsum[NR] = 0; sum += $3
      if (NR == 1) root = $2
      while (top > 0 && sdep[top] >= dep) top--
      if (top > 0) childsum[stk[top]] += $2
      stk[++top] = NR; sdep[top] = dep }
    END { for (i = 1; i <= NR; i++) if (clade[i] != taxon[i] + childsum[i]) bad++
          if (bad || sum != root) { printf "rows not adding up %d, taxon sum %d, root %d\n", bad, sum, root; exit 1 } }' "$1"
}

echo "== kraken-report read counts"
# Codex's counterexample: reads shared by three sibling strains (genome lengths 1M, 1M and 2M).
cp "$X/nodes.dmp" "$T/nodes3.dmp"; cp "$X/names.dmp" "$T/names3.dmp"
printf '900001\t|\t446\t|\tstrain\t|\n' >> "$T/nodes3.dmp"
printf '900001\t|\tStrain_C\t|\t\t|\tscientific name\t|\n' >> "$T/names3.dmp"
printf '272624 1000000\n297246 1000000\n900001 2000000\n' > "$T/sizes3.tsv"
for n in 1 3; do
  {
    printf 'readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n'
    for r in $(seq 1 $n); do
      for t in 272624 297246 900001; do printf 'r%d\ts\t%d\t100\t100\t100\t100\t3\n' "$r" "$t"; done
    done
  } > "$T/shared$n.tsv"
done
for em in "" "--em-reads"; do
  name=${em:-default}
  "$Q" -x "$T/ex" -c "$T/class.tsv" --output-format 3 $em > "$T/ex.$name.f3.tsv" 2>/dev/null
  if msg=$(kreport_adds_up "$T/ex.$name.f3.tsv"); then ok "example, $name: counts add up"; else bad "example, $name: $msg"; fi
  for n in 1 3; do
    "$Q" --taxonomy-tree "$T/nodes3.dmp" --name-table "$T/names3.dmp" --size-table "$T/sizes3.tsv" \
      -c "$T/shared$n.tsv" --output-format 3 $em > "$T/shared$n.$name.f3.tsv" 2>/dev/null
    root=$(awk -F'\t' 'NR==1{print $2}' "$T/shared$n.$name.f3.tsv")
    if msg=$(kreport_adds_up "$T/shared$n.$name.f3.tsv") && [ "$root" = "$n" ]; then
      ok "$n read(s) shared by three strains, $name: counts add up, root $n"
    else bad "$n read(s) shared by three strains, $name: ${msg:-root $root}"; fi
  done
done
# Without --em-reads the rows and the percentage column are the same as before; only read counts are made whole.
"$Q" -x "$T/ex" -c "$T/class.tsv" > "$T/ex.f0.tsv" 2>/dev/null
rows=$(awk -F'\t' 'NR>1 && $5>=1e-6' "$T/ex.f0.tsv" | wc -l)
[ "$(wc -l < "$T/ex.default.f3.tsv")" -ge "$rows" ] && ok "default kraken-report lists every taxon with reads" \
  || bad "default kraken-report lists every taxon with reads"

echo "== protein index read weights"
"$Q" -x "$T/ex" -c "$T/class.tsv" > "$T/nt.tsv" 2>/dev/null
"$Q" -x "$T/aa" -c "$T/class_aa.tsv" > "$T/aa_x.tsv" 2> "$T/aa_x.log"
"$Q" --taxonomy-tree "$X/nodes.dmp" --name-table "$X/names.dmp" -c "$T/class_aa.tsv" --protein > "$T/aa_flag.tsv" 2>/dev/null
"$Q" --taxonomy-tree "$X/nodes.dmp" --name-table "$X/names.dmp" -c "$T/class_aa.tsv" > /dev/null 2> "$T/aa_noflag.log"
"$Q" -x "$T/ex" -c "$T/class.tsv" > /dev/null 2> "$T/nt.log"
# emNumReadsWeighted at the root: with amino-acid hits it should be about what nucleotide hits give, not ~0.
w() { awk -F'\t' 'NR>1 && $2==1 {print $9}' "$1"; }
wnt=$(w "$T/nt.tsv"); wx=$(w "$T/aa_x.tsv"); wf=$(w "$T/aa_flag.tsv")
awk -v a="$wnt" -v b="$wx" -v c="$wf" 'BEGIN{exit !(a > 1 && b > a / 2 && b < a * 2 && c > a / 2 && c < a * 2)}' \
  && ok "root emNumReadsWeighted: nucleotide $wnt, -x protein $wx, --protein $wf" \
  || bad "root emNumReadsWeighted: nucleotide $wnt, -x protein $wx, --protein $wf"
grep -q "Protein index" "$T/aa_x.log" && ok "protein index detected from .4.cfr" || bad "protein index detected from .4.cfr"
grep -q "Warning: no hit covers more than a third" "$T/aa_noflag.log" && ok "warns when amino-acid hits are read as nucleotide" \
  || bad "warns when amino-acid hits are read as nucleotide"
grep -q "Warning: no hit covers" "$T/nt.log" && bad "no warning for a nucleotide index" || ok "no warning for a nucleotide index"

if [ $fail -eq 0 ]; then echo "all tests passed"; else echo "some tests FAILED"; fi
exit $fail
