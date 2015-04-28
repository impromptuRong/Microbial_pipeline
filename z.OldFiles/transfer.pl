#! usr/bin/perl
#! use strict; use warnings;
###########  Read in List files  ##############
$file = $ARGV[0];
$cutoff = $ARGV[1];

my %Bac;
open(RRDP,"<$file.shhh.trim.rdp.output") or die "Cannot open old $file.rdp.output\n";
open(ROUT,">$file.shhh.trim.rdp.taxonomy.$cutoff") or die "Cannot open old $file.rdp.rdp.taxonomy\n";
local $/ = "\n";
foreach my $line(<RRDP>){
	chomp($line);
	$line =~ s/ /_/g;
	my @info = split(/\t/,$line);
	my $seqname = $info[0];

	my @index = grep {$info[$_] eq "domain"} 0..$#info;
	my $domain = @index?$info[$index[0]-1]:"unclassified";
	my $dscore = @index?$info[$index[0]+1]:0;
	$dscore = 100*$dscore;
	
	my @index = grep {$info[$_] eq "phylum"} 0..$#info;
	my $phylum = @index?$info[$index[0]-1]:"unclassified";
	my $pscore = @index?$info[$index[0]+1]:0;
	$pscore = 100*$pscore;

	my @index = grep {$info[$_] eq "class"} 0..$#info;
	my $class = @index?$info[$index[0]-1]:"unclassified";
	my $cscore = @index?$info[$index[0]+1]:0;
	$cscore = 100*$cscore;

	my @index = grep {$info[$_] eq "order"} 0..$#info;
	my $order = @index?$info[$index[0]-1]:"unclassified";
	my $oscore = @index?$info[$index[0]+1]:0;
	$oscore = 100*$oscore;

	my @index = grep {$info[$_] eq "family"} 0..$#info;
	my $family = @index?$info[$index[0]-1]:"unclassified";
	my $fscore = @index?$info[$index[0]+1]:0;
	$fscore = 100*$fscore;

	my @index = grep {$info[$_] eq "genus"} 0..$#info;
	my $genus = @index?$info[$index[0]-1]:"unclassified";
	my $gscore = @index?$info[$index[0]+1]:0;
	$gscore = 100*$gscore;

	print ROUT $seqname,"\t";
	print ROUT $dscore>$cutoff?"$domain($dscore);":"unclassified;";
	print ROUT $pscore>$cutoff?"$phylum($pscore);":"unclassified;";
	print ROUT $cscore>$cutoff?"$class($cscore);":"unclassified;";
	print ROUT $oscore>$cutoff?"$order($oscore);":"unclassified;";
	print ROUT $fscore>$cutoff?"$family($fscore);":"unclassified;";
	print ROUT $gscore>$cutoff?"$genus($gscore);":"unclassified;";
	print ROUT "\n";
}
close RRDP;
close ROUT;

