#! usr/bin/perl
#! use strict; use warnings;

$filename = $ARGV[0];

open(FILE, "./$filename") or die $!;
local $/ = "\n\n\n\n";
my $line = <FILE>;
#$line =~ s/(\n|\t|\s)//g;
close FILE;

my @levelarray=("genus","family","order","class","phylum");
my $info="";
foreach my $level(@levelarray){
	if($filename =~ /$level/){$info=$level;}
}
print $info,"\n";
open(DIC, "./dic.$info.csv") or die "Failed to open dic.$info.csv\n";
my $dic = <DIC>;
my @trans = split(/\n/,$dic);
for(my $i=scalar(@trans);$i>1;$i--){
	my @tmp = split(/,/, $trans[$i-1]);
	my $pattern = $tmp[1];
	my $instead = $tmp[0];
	print $pattern,"\t",$instead,"\n";
#	$line =~ s/\Q$pattern\E/$instead/g;
	$line =~ s/$pattern/$instead/g;
}
close DIC;

open(OUT, ">./$filename.trans");
print OUT $line;
close OUT;

