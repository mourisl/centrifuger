#!/bin/perl

use strict ;
use warnings ;

use Getopt::Long ;

my $usage = "perl ./download-gtdb.pl [options]:\n".
  "options:\n".
  "\t-o STR: output prefix [gtdb]\n".
  "\t-t: number of threads [1]\n".
  "\t--names STR: the NCBI names.dmp file\n".
  "\t-h: print this message and exit\n".
  "\n" ;


sub PrintLog
{
  print STDERR "[".localtime()."] ".$_[0]."\n" ;
}

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my $outputPrefix = "gtdb" ;
my $printHelpAndDie = 0 ;
my $threadCnt = 1 ;
my $ncbiNameDmp = "" ;

GetOptions
(
  "o=s" => \$outputPrefix,
  "t=i" => \$threadCnt,
  "h" => \$printHelpAndDie,
  "names=s"=>\$ncbiNameDmp 
) ;

if ($printHelpAndDie)
{
  die "$usage" ; 
}

my $ftp="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/" ;
my $genomeTarFile="gtdb_genomes_reps.tar.gz" ;
system_call("curl -o $genomeTarFile $ftp/genomic_files_reps/$genomeTarFile") ;
system_call("tar -xzf $genomeTarFile") ;

# Get the version
system_call("curl -o gtdb_version.txt $ftp/VERSION.txt") ;
open FP, "gtdb_version.txt" ;
my $line = <FP> ;
chomp $line ;
my $gtdbVersion = substr($line, 1) ;
close FP ;

# Create the metadata file that combines bac and archaea
system_call("curl $ftp/bac120_metadata.tsv.gz | gzip -cd > ${outputPrefix}_meta.tsv") ;
system_call("curl $ftp/ar53_metadata.tsv.gz | gzip -cd | grep -v \"^accession\" >> ${outputPrefix}_meta.tsv") ;

my $createDmpOptions = " -d gtdb_genomes_reps_r$gtdbVersion -m ${outputPrefix}_meta.tsv -o $outputPrefix -t $threadCnt" ;
$createDmpOptions .= " --names $ncibNameDmp" if ($ncibNameDmp ne "") ;

system_call("perl gtdb-create-dmp.pl $createDmpOptions") ;
