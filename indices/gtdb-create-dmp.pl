#!/bin/perl

use strict ;
use warnings ;

use Cwd qw(abs_path) ;
use Getopt::Long ;
use threads ;
use threads::shared ;

# Create the dump files that are necessary for centrifuger-build from GTDB files
sub AccessionToSubdir
{
  my $accession = $_[0] ;
  return substr($accession, 0, 3)."/".substr($accession, 4, 3).
        "/".substr($accession, 7, 3)."/".substr($accession, 10, 3) ;
}

sub GetGenomeFilePath
{
  my $dir = $_[0] ; 
  my $accession = $_[1] ; #GCF_000657795.2
  my $subdir = AccessionToSubdir($accession)  ;

  return $dir."/database/".$subdir."/".$accession."_genomic.fna.gz" ;
}

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

my $usage = "Usage: gtdb-create-dmp.pl [OPTIONS]\n".
  "\t-d STR: directory of GTDB untarred representative sequence\n".
  "\t-m STR: GTDB metadata file\n".
  "\t-o STR: output prefix [gtdb_]\n".
  #"\t-t INT: number of threads to use [1]\n".
  "\t--names STR: NCBI's names.dmp file. If not given, using non-NCBI taxid to represent intermediate nodes.\n".
  "\t--taxIDStart INT: the tax ID number start to count. [10000000]\n".
  "\t--generateSeqId2TaxId: generate seqid_to_taxid.map conversion file.\n"
  ;

die "$usage\n" if (@ARGV == 0) ;

# Generate names.dmp and nodes.dmp files using metadata file
PrintLog("Generate the dmp files for nodes and names, and the mapping from file path to tax ID.") ;

my $ncbiNodeDmp = "" ;
my $ncbiNameDmp = "" ;
my $outputPrefix = "gtdb" ;
my $genomeDir = "" ;
my $metaFile = "" ;
my $novelTaxId = 10000000 ;
my $numThreads = 1 ;
my $generateSeqIdMap = 0 ;

GetOptions(
  "o=s" => \$outputPrefix,
  "d=s" => \$genomeDir,
  "m=s" => \$metaFile,
  #"t=i" => \$numThreads, 
  #"nodes=s" => \$ncbiNodeDmp,
  "names=s" => \$ncbiNameDmp,
  "taxIDStart=i"=> \$novelTaxId,
  "generateSeqId2TaxId" => \$generateSeqIdMap 
) ;
my $fullGenomeDir = abs_path($genomeDir) ;

my $i ;
my @cols ;

my %taxRankCodeToFull = ("d"=>"domain", "p"=>"phylum", "c"=>"class",
         "o"=>"order", "f"=>"family", "g"=>"genus", "s"=>"species", "x"=>"no rank" ) ;

# Collect the NCBI names information
#3041336	|	Stipitochalara	|		|	scientific name	|
my %ncbiNamesToTaxId ;
if ($ncbiNameDmp ne "")
{
  open FP, $ncbiNameDmp ;
  while (<FP>)
  {
    chomp ;
    @cols = split /\t/ ;
    next if ($cols[6] ne "scientific name") ;

    my $name = $cols[2] ;
    #$name =~ s/\s/_/g ; 
    $ncbiNamesToTaxId{$name} = $cols[0] ;
  }
  close FP ;
}


open FPoutNames, ">${outputPrefix}_names.dmp" ;
open FPoutNodes, ">${outputPrefix}_nodes.dmp" ;
open FPoutFileToTaxid, ">${outputPrefix}_fname_to_taxid.map" ;
open FPoutFileList, ">${outputPrefix}_file.list" ;
open FPmeta, $metaFile ;

print FPoutNodes "1\t|\t1\t|\tno rank\t|\n" ;
print FPoutNames "1\t|\troot\t|\tscientific name\t|\n" ;

my $header = <FPmeta> ;
my %colNames ; 

# We use GTDB meta instead of GTDB taxonomy is the meta directly provides the 
#   representative genomes. So we can direclty generate the file list.
@cols = split /\t/, $header ; 
for ($i = 0 ; $i < scalar(@cols) ; ++$i)
{
  $colNames{ $cols[$i] } = $i ;
}

my %fileNameToTaxId ;
my %accessionToTaxId ;
my %nodesToPrint ;
my %taxIdRank ;
my %namesToPrint ;
my %newNamesToTaxId ; # the names that is not in NCBI name table. 

while (<FPmeta>)
{
  @cols = split /\t/ ;
  next if ($cols[ $colNames{"gtdb_representative"}] ne "t") ;
  
  my $accession = substr($cols[ $colNames{"accession"} ], 3) ;

  #if (! -e GetGenomeFilePath($genomeDir, $accession) )
  #{
  #  print("Warning: $accession not exists\n") ;
  #  next ;
  #}

  # Use ncbi_taxid may create issues as GTDB may reassign 
  #   a strain to another species, so the ncib_taxid still point
  #   towards the original species/strain, but the
  #   lineage string is totally different.
  #my $taxid = $cols[ $colNames{"ncbi_taxid"} ] ; 
  my $taxid = 1 ;
  my $lineage = $cols[ $colNames{"gtdb_taxonomy"} ] ;
  
  #next if ($taxid eq "none") ;

  #d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Bordetella;s__Bordetella pseudohinzii
  my @lineageFields = split /;/, $lineage ;

  my $j ;
  my $parentTid = 1 ;
  my $ltid = 1 ;
  for ($j = 0 ; $j < scalar(@lineageFields) ; ++$j)
  {
    my @cols2 = split /__/, $lineageFields[$j] ;

    if (defined $ncbiNamesToTaxId{$cols2[1]})
    {
      $ltid = $ncbiNamesToTaxId{$cols2[1]} ;
    }
    elsif (defined $newNamesToTaxId{$lineageFields[$j]}) # Sometimes GTDB put the same name across multiple lineage trees, so we need to add the prefix like o__ to distinguish them
    {
      $ltid = $newNamesToTaxId{$lineageFields[$j]} ;
    }
    else
    {
      $ltid = $novelTaxId ;
      #$ncbiNamesToTaxId{$lineageFields[$j]} = $ltid ;
      $newNamesToTaxId{$lineageFields[$j]} = $ltid ; 
      ++$novelTaxId ;
    }

    $taxid = $ltid if ($j == scalar(@lineageFields) - 1) ;

    if (defined $nodesToPrint{$ltid} && $nodesToPrint{$ltid} != $parentTid)
    {
      die "A conflict of lineage information is found when processing $lineage\n" ;
    }

    $nodesToPrint{$ltid} = $parentTid ;
    $taxIdRank{$ltid} = $cols2[0] ;
    $namesToPrint{$ltid} = $cols2[1]  ;

    $parentTid = $ltid ;
  }

  $accessionToTaxId{$accession} = $taxid ;
  my $fileName = GetGenomeFilePath($fullGenomeDir, $accession) ;
  $fileNameToTaxId{$fileName} = $taxid ;
  print FPoutFileToTaxid "$fileName\t$taxid\n" ;
  print FPoutFileList "$fileName\n" ;
}

foreach my $tid (keys %nodesToPrint)
{
  print FPoutNodes "$tid\t|\t".$nodesToPrint{$tid}."\t|\t".$taxRankCodeToFull{$taxIdRank{$tid}}."\t|\n" ;
  print FPoutNames "$tid\t|\t".$namesToPrint{$tid}."\t|\tscientific name\t|\n" ;
}

close FPmeta ;
close FPoutNames ;
close FPoutNodes ;
close FPoutFileToTaxid ;
close FPoutFileList ;

# Iterate through the genome files to generate the seqid map file
if ($generateSeqIdMap == 0)
{
  exit(0) ;
}

PrintLog("Generate the seq ID to tax ID mapping file.") ;
#cat gtdb_file.list | xargs -I {} sh -c 'gzip -cd {} | grep "^>" | while read header; do echo "$header {}"; done' 
system_call("cat ${outputPrefix}_file.list | xargs -I {} sh -c 'gzip -cd {} | grep \"^>\" | while read header; do echo \"\$header {}\"; done' > ${outputPrefix}_faheader_to_file.out") ;

my %seqIdMap ;
open FP, "${outputPrefix}_faheader_to_file.out" ;
while (<FP>)
{
  chomp ;
  my @cols = split /\s+/, $_ ;
  my $colNum = scalar(@cols) ;
  
  my $seqId = substr($cols[0], 1) ;
  my $fileName = $cols[ $colNum - 1] ;
  my $taxId = $fileNameToTaxId{$fileName} ;

  $seqIdMap{$seqId} = $taxId ;
}
close FP ;
unlink("${outputPrefix}_faheader_to_file.out") ;

open FPout, ">${outputPrefix}_seqid_to_taxid.map" ;
foreach my $seqId (keys %seqIdMap)
{
  print FPout "$seqId\t".$seqIdMap{$seqId}."\n" ;
}
close FPout ;

PrintLog("Done.") ;
