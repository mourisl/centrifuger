#!/bin/perl

use strict ;
use warnings ;

use Getopt::Long ;

my $usage = "perl ./gtdb-download.pl [options]:\n".
  "options:\n".
  "\t-o STR: output prefix [gtdb]\n".
  "\t--names STR: the NCBI names.dmp file\n".
  "\t--generateSeqId2TaxId: generate seqid_to_taxid.map conversion file.\n".
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
my $ncbiNameDmp = "" ;
my $generateSeqIdMap = 0 ;

GetOptions
(
  "o=s" => \$outputPrefix,
  "h" => \$printHelpAndDie,
	"generateSeqId2TaxId" => \$generateSeqIdMap,
  "names=s"=>\$ncbiNameDmp
) ;

if ($printHelpAndDie)
{
  die "$usage" ; 
}

my $ftp="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/" ;
my $genomeTarFile="${outputPrefix}_genomes_reps.tar.gz" ;
system_call("curl -o $genomeTarFile $ftp/genomic_files_reps/$genomeTarFile") ;
system_call("tar -xzf $genomeTarFile") ;

# Get the version
system_call("curl -o ${outputPrefix}_version.txt $ftp/VERSION.txt") ;
open FP, "${outputPrefix}_version.txt" ;
my $line = <FP> ;
chomp $line ;
my $gtdbVersion = substr($line, 1) ;
close FP ;

# Create the metadata file that combines bac and archaea
system_call("curl $ftp/bac120_metadata.tsv.gz | gzip -cd > ${outputPrefix}_meta.tsv") ;
system_call("curl $ftp/ar53_metadata.tsv.gz | gzip -cd | grep -v \"^accession\" >> ${outputPrefix}_meta.tsv") ;

my $createDmpOptions = " -d gtdb_genomes_reps_r$gtdbVersion -m ${outputPrefix}_meta.tsv -o $outputPrefix" ;
$createDmpOptions .= " --names $ncbiNameDmp" if ($ncbiNameDmp ne "") ;
$createDmpOptions .= " --generateSeqId2TaxId" if ($generateSeqIdMap) ;

system_call("perl gtdb-create-dmp.pl $createDmpOptions") ;
