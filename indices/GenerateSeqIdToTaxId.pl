#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: perl GenerateSeqIdToTaxId.pl genomes.fa.gz accession2taxid.gz > seqid_to_taxid.map\n" if (@ARGV == 0) ;

open FP, "gzip -cd $ARGV[0] |" ;
my %seqids ;
while (<FP>)
{
  next if (!/^>/) ;
  $seqids{ substr((split)[0], 1)} = 1 ;
}
close FP ;

open FP, "gzip -cd $ARGV[1] |" ;
<FP> ; # header
while (<FP>)
{
  chomp ;
  my @cols = split ;
  print $cols[1]."\t".$cols[2]."\n" if (defined $seqids{$cols[1]}) ;
}
close FP ;
