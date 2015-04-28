#! usr/bin/perl
#! use strict; use warnings;
use POSIX qw/ceil/;
use POSIX qw/floor/;
#! use Math::Round qw/round/;
###########  Read in List files  ##############
$file = $ARGV[0];

open(RDP,"<$file") or die "Cannot open old $file\n";
my @tmp = split(/\./,$file);
my $cutoff = $tmp[-1];
my $length = ceil((100-$cutoff)/10)+3;
my @dcount = ((0)x$length);
local $/ = "\n";
foreach my $line(<RDP>){
	chomp($line);
	$line =~ s/ /_/g;
	my ($seqname,$info) = split(/\t/,$line);
	my @tmp = split(/;/,$info);
	my $dinfo = $tmp[0];

	if($dinfo eq "unclassified"){$dcount[-1]++;}
	else{
		my @tmp = split(/\(|\)/,$dinfo);
		my $domain = $tmp[0];
		my $dscore = $tmp[1];
		if($domain eq "Bacteria"){
			my $i = ceil((100-$dscore)/10);
			$dcount[$i]++;
		}
		else{
			$dcount[-2]++;
		}
	}
}

for(my $i=0;$i<($length-3);$i++){
	my $s = 100-10*$i;
	print $s,"\t\t",$dcount[$i],"\n";
}
print $cutoff, "\t\t", $dcount[-3],"\n";
print "unclassified\t", $dcount[-1],"\n";
print "non-bac\t\t", $dcount[-2],"\n";

close RDP;

