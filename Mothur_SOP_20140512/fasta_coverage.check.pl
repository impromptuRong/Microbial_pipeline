#! usr/bin/perl
#! use strict; use warnings;

$query = $ARGV[0];
$template = $ARGV[1];

open(TFSA, "<./$template") or die "Failed to open $template\n";
open(QFSA, "<./$query") or die "Failed to open $query\n";
open(OUT, ">./coverage.check.csv") or die "Failed to open Output file\n";

############    Read Template Fasta File    ############
my %tfsa;
local $/ = ">";
$_ = <TFSA>;
while(<TFSA>){
	chomp;
	my @array = split(/\n+/);
	my ($seqname,$junk) = split(/\s+/,shift(@array));
	my $seq = join("",@array);
	$seq =~ s/\s|-|\.//g;
	$seq =~ s/U/T/g;
	$seq =~ s/u/t/g;
	$tfsa{$seqname} = $seq;
}
close TFSA;

############    Read Query Fasta File    ############
local $/ = ">";
$_ = <QFSA>;
while(<QFSA>){
	chomp;
	my @array = split(/\n+/);
	my ($seqname,$junk) = split(/\s+/,shift(@array));
	my $seq = join("",@array);
	$seq=~ s/\s|-|\.//g;
	$seq =~ s/U/T/g;
	$seq =~ s/u/t/g;
	print OUT $seqname,",",length($seq),",",length($tfsa{$seqname}),",",length($seq)/length($tfsa{$seqname}),"\n";
}
close QFSA;
close OUT;
