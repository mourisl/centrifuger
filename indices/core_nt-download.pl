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
  "\t--clean: whether remove the blastdb files [keep]\n"
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
my $clean = 0 ;

GetOptions
(
  "o=s" => \$outputPrefix,
  "blast=s" => \$blastPath
  "clean" => \$clean
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

system_call("${blastPath}/blastdbcmd -entry all -db core_nt -out /dev/stdout | gzip -c > ${outputPrefix}.fa.gz") ;

system_call("rm core_nt.*.tar.gz") ;
if ($clean == 1)
{
  if ($outputPrefix ne "core_nt")
  {
    system_call("rm core_nt.*") ;
  }
  else
  {
    system_call("ls core_nt.* | grep -v core_nt.fa.gz | xargs -I{} rm {}") ;
  }
}

# User need to download: wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# and run centrifuger-download taxonomy to get the nodes.dmp and names.dmp file



