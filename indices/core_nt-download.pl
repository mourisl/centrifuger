#!/usr/bin/env perl

use warnings ;
use strict ;

use Cwd qw(abs_path) ;
use Getopt::Long ;
use threads ;
use threads::shared ;

my $usage = "perl ./core_nt-download.pl [options]:\n".
  "\t--blast STR: path to blast bin folder hold blastdbcmd []\n". 
  "\t--accession2taxid STR: the accesion2taxid file\n".
  "\t-o STR: output prefix [core_nt]\n".
  "\t--noclean: whether remove the blastdb files and intermediate files [remove]\n"
  if (@ARGV == 0) ;

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

my $outputPrefix = "core_nt" ;
my $blastPath = "" ;
my $noclean = 0 ;

GetOptions
(
  "o=s" => \$outputPrefix,
  "blast=s" => \$blastPath
  "noclean" => \$noclean
) ;

# get the download list
system_call("wget https://ftp.ncbi.nih.gov/blast/db/core_nt-nucl-metadata.json") ;

open FP, "core_nt-nucl-metadata.json" ;
my @tars ;
while (<FP>)
{
  chomp ;
  my $line = $_ ;

  if ($line =~ /(ftp:\/\/ftp.ncbi.nlm.nih.gov\/blast\/db\/core_nt.[0-9]+.tar.gz)/)
  {
    print $1, "\n" ;  
    push @tars, $1 ;
  }
  
}
close FP ;

# download each sample
foreach my $link (@tars)
{
  system_call("wget $link") ;
}

# uncompress
foreach my $link (@tars)
{
  if ($link =~ /(core_nt.[0-9]+.tar.gz)/) 
  {
    system_call("tar -xzf $1") ;   
  }
  else
  {
    die "Wrong link format $link\n" ;
  }
}

# Some how this gzipped output always fail on my system.
# It might be fine, as the dustmasker requires uncomrpessed input.
#system_call("${blastPath}/blastdbcmd -entry all -db core_nt -out - | gzip -c > ${outputPrefix}.fa.gz") ;
system_call("${blastPath}/blastdbcmd -entry all -db core_nt -out ${outputPrefix}.fa") ;

system_call("rm core_nt.*.tar.gz") ;
if ($noclean == 0)
{
  if ($outputPrefix ne "core_nt")
  {
    system_call("rm core_nt.*") ;
  }
  else
  {
    system_call("ls core_nt.* | grep -v core_nt.fa. | xargs -I{} rm {}") ;
  }
}

# User need to download: wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# and run centrifuger-download taxonomy to get the nodes.dmp and names.dmp file
system_call("${blastPath}/dustmasker -in ${outputPrefix}.fa -outfmt fasta | gzip -c > ${outputPrefix}_dustmasker.fa.gz") ;

if ($noclean == 0)
{
  system_call("rm ${outputPrefix}.fa") ;
}
# Lastly to run "perl GenerateSeqIdToTaxId.pl core_nt_dustmasker.fa.gz nucl_gb.accession2taxid.gz > core_nt_seqid_to_taxid.map"

