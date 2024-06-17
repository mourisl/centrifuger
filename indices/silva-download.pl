#!/usr/bin/env perl

use strict ;
use warnings ;

my $help = "usage: download-silva.pl [OPTIONS]:\n".
  "\t-v STRING: SILVA_version [138.1]\n".
  "\t-o STRING: output_folder [./]\n".
  "\t--subunit STRING: use SSU or LSU [SSU]\n".
  "\t--NR99 INT: use NR99 genome or not: 0 for not, 1 for yes [1]\n".
  "\t-h: print help message and quit.\n";

my $i ;
my $silvaVer = "138.1" ;
my $outputDir = "./" ;
my $useNR99 = 1 ;
my $subunit = "SSU" ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
}

for ($i = 0 ; $i < scalar(@ARGV) ; ++$i)
{
  if ($ARGV[$i] eq "-v")
  {
    $silvaVer = $ARGV[$i + 1] ;
    ++$i ;
  }
  elsif ($ARGV[$i] eq "-o")
  {
    $outputDir = $ARGV[$i + 1] ;
    ++$i ;
  }
  elsif ($ARGV[$i] eq "--subunit")
  {
    $subunit = $ARGV[$i + 1] ;
    ++$i ;
  }
  elsif ($ARGV[$i] eq "--NR99")
  {
    $useNR99 = $ARGV[$i + 1] ;
    ++$i ;
  }
  elsif ($ARGV[$i] eq "-h")
  {
    die "$help" ;
  }
}

my $underscoreVer = $silvaVer ;
$underscoreVer =~ s/\./_/g ;

#https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz
my $weblink = "https://www.arb-silva.de/fileadmin/silva_databases/release_$underscoreVer/Exports" ;
my $prefix = "" ;

# download and create the taxonomy tree and name file
$prefix = "tax_slv_".lc($subunit)."_$silvaVer" ;
system_call("wget $weblink/taxonomy/${prefix}.txt.gz") ;

# Archaea;Aenigmarchaeota;  11084 phylum    123
open FP, "zcat ${prefix}.txt.gz |" ;
my %nameMap ; 
open FPnames, ">$outputDir/names.dmp" ;
# Get the names
print FPnames "1\t|\troot\t|\tscientific name\t|\n" ;
while (<FP>)
{
  chomp ;
  my @cols = split /\t/, $_ ;
  my @nameCols = split /;/, $cols[0] ;
  my $name = $nameCols[scalar(@nameCols) - 1] ;
  my $tax = $cols[1] ;
  $nameMap{$name} = $tax ;

  print FPnames "$tax\t|\t$name\t|\tscientific name\t|\n" ;
}
close FPnames ;
close FP ;

open FP, "zcat ${prefix}.txt.gz |" ;
open FPnodes, ">$outputDir/nodes.dmp" ;
print FPnodes "1\t|\t1\t|\tno rank\t|\n" ;
while (<FP>)
{
  chomp ;
  my @cols = split /\t/, $_ ;
  my @nameCols = split /;/, $cols[0] ;
  my $tax = $cols[1] ;
  my $parentName ;
  $parentName = $nameCols[scalar(@nameCols) - 2] if (scalar(@nameCols) > 1) ;
  my $parent ;
  $parent = $nameMap{$parentName} if (defined $parentName && defined $nameMap{$parentName}) ;
  $parent = 1 if (!defined $parent) ;
  print FPnodes "$tax\t|\t$parent\t|\t".$cols[2]."\t|\n" ;
}
close FPnodes ;
close FP ;

unlink "${prefix}.txt.gz"  ;

# download the seqid_to_taxid map file
system_call("wget $weblink/taxonomy/${prefix}.acc_taxid.gz") ;
system_call("zcat ${prefix}.acc_taxid.gz > $outputDir/silva_seqid_to_taxid.map") ;
unlink "${prefix}.acc_taxid.gz" ;

# download the genome file
$prefix = "SILVA_".$silvaVer."_".$subunit."Ref_" ;
if ($useNR99 == 1)
{
  $prefix .= "NR99_" ;
}

system_call("wget ${weblink}/${prefix}tax_silva.fasta.gz") ;
open FP, "zcat ${prefix}tax_silva.fasta.gz |" ;
open FPout, "| gzip -c > $outputDir/silva_seq.fa.gz" ;
my %seqToTax ;
while (<FP>)
{
  my $line = $_ ;
  chomp $line ;
  if ($line =~ /^>/)
  {
    my @cols = split /\s/, $line ;
    my $seqId = $cols[0] ;
    print FPout "$seqId\n" ;
  }
  else # change the nucleotide sequence from U to T
  {
    $line =~ s/U/T/g ;
    print FPout "$line\n" ;
  }
}
close FP ;
close FPout ;
unlink "${prefix}tax_silva.fasta.gz" ;

