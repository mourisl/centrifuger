#!/bin/perl

use strict ;
use warnings ;

# seqid not in accession file will be assigned to taxonomy ID 1.
die "perl a.pl seqid.list accession_taxid.map > seqid_to_taxid.map\n" if (@ARGV == 0) ;

my %accessionMap ;
open FP, $ARGV[1] ;
while (<FP>)
{
  chomp ;
  my @cols = split ;
  $accessionMap{$cols[0]} = $cols[2] ; # No harm to include the header
}
close FP ;

open FP, $ARGV[0] ;
while (<FP>)
{
  chomp ;
  my $line = $_ ;
  my @cols = split /\./, $line ;
  if (defined $accessionMap{$cols[0]})
  {
    print("$line\t", $accessionMap{$cols[0]}, "\n") ;
  }
  else
  {
    print("$line\t1\n") ;
  }
}
close FP ;
