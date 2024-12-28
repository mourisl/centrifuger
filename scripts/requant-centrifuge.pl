#!/usr/bin/env perl

use strict ;
use warnings ;

die "Usage: $0 centrifuger_path centrifuge_path centrifuge_index_prefix centrifuge_class_out > new_report.tsv \n" if (@ARGV == 0) ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
}

my $prefix = "tmp" ;
my $nodes = "${prefix}_nodes.out" ;
my $names = "${prefix}_names.out" ;
my $sizes = "${prefix}_sizes.out" ;

my $cfInspect = $ARGV[0]."/centrifuge-inspect" ;
my $cfIndex = $ARGV[2] ;

system_call("$cfInspect --taxonomy-tree $cfIndex > $nodes") ;
system_call("$cfInspect --name-table $cfIndex | awk '{print \$1\"\\t|\\t\"\$2\"\\t|\\tscientific name\"}'> $names") ;
system_call("$cfInspect --size-table $cfIndex > $sizes") ;

my $cfrQuant = $ARGV[1]."/centrifuger-quant" ;
my $classification = $ARGV[3] ;

system_call("$cfrQuant --taxonomy-tree $nodes --name-table $names --size-table $sizes -c $classification") ;

unlink $nodes ;
unlink $names ;
unlink $sizes ;
