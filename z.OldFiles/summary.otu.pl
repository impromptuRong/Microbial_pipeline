#! usr/bin/perl
#! use strict; use warnings;
use List::Util qw/shuffle/;

$file = $ARGV[0];
my %S;
open(GRP, "<./$file.final.groups") or die "Failed to open Groups File\n";
local $/ = "\n";
my %barhash;
while(<GRP>){
	chomp;
	my ($seqname,$bar) = split(/\t/);
	$S{$seqname}{"barcode"} = $bar;
	$barhash{$bar} = 1;		
}
close GRP;
my @bararray = sort keys %barhash;

open(TAX, "<./$file.final.unique.an.0.03.cons.taxonomy") or die "Failed to open Taxonomy File\n";
my %info;
$_ = <TAX>;
foreach (<TAX>){
	chomp;
	my @tmp = split(/\t/);
	my $otu = sprintf("Otu%04d",$tmp[0]);
	my $total = $tmp[1];
	$tmp[2] =~ s/;/,/g;
	my $index = join("",$tmp[2],$total);
	$info{$otu} = $index;
}

open(OUT, ">./summary.otu.csv") or die "Failed to open Taxonomy File\n";
print OUT "Group,Domain,Phylum,Class,Order,Family,Genus,Barcode"; 
foreach my $barcode(@bararray){
	print OUT ",",$barcode;
} 

open(OTU, "<./$file.final.unique.an.0.03.otu") or die "Failed to open OTUlist File\n";

local $/ = "\n";
while(<OTU>){
	chomp;
	my ($otu,$tmp) = split(/\t/);
	my $otu = sprintf("Otu%04d",$otu);
	print OUT "\n",$otu,",",$info{$otu};
	my @seqarray = split(/,/,$tmp);
	my %otuline;
	foreach my $seqname(@seqarray){
		$otuline{$S{$seqname}{"barcode"}} += 1;
	}
	foreach my $barcode(@bararray){
		print OUT ",";
		exists($otuline{$barcode})? print OUT $otuline{$barcode} : print OUT 0;
	}
}
close OTU;
close TAX;
close OUT;

