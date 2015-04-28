#! usr/bin/perl
#! use strict; use warnings;
use List::Util qw/shuffle/;

########## perl summary.genera.pl 
if(scalar @ARGV==2){
	$depth=$ARGV[0];
	$permu=$ARGV[1];
}
elsif(scalar @ARGV==1){
	$depth=$ARGV[0];
	$permu=1;
#	$cutoff=if($ARGV[1]>1)?$ARGV[1]:($ARGV[1]*100);
}

else{
	$depth="all";
	$permu=1;
}
print $depth,"\t",$permu,"\n";

system("rm -rf $depth.$permu");
system("mkdir $depth.$permu");
###########  Read in List files  ##############
opendir(SEP,"./barcode.shhh.fasta") or die "Cannot Open fasta Dir!";
my @filefolder = readdir(SEP);

my %B;
open(AAH, ">genera");
foreach my $file (@filefolder){
if($file =~ /fasta/){
	open(OFSA,"<./barcode.shhh.fasta/$file") or die "Cannot open shhh $file\n";
	my @tmp = split(/\./,$file);
	pop(@tmp);
	my $bar = join(".",@tmp);
	local $/ = ">";
	$_ = <OFSA>;
	while(<OFSA>){
		chomp;
		my @array = split(/\n+/);
		my @tmp = split(/\t+/,shift(@array));
		my $seqname = $tmp[0];
		my $genera = ($tmp[-1]eq"unclassified")?"unclassified": $tmp[-2];
		my $seq = join("",@array);
		$seq=~ s/\s|-|\.//g;
		$B{$bar}{$seqname}{"seq"} = $seq;
		$B{$bar}{$seqname}{"info"} = $genera;
print AAH $bar,"\t",$seqname,"\t",$genera,"\n";
	}
	close OFSA;
}
}
close AAH;
for(my $n=1;$n<=$permu;$n++){
        system("rm -rf ./$depth.$permu/$n");
        system("mkdir ./$depth.$permu/$n");
print "Selection: $n\n";
	my %P;
	foreach $bar (keys %B){
		my @total = keys %{$B{$bar}};
		my @select;
###########  Random Selection  ##############
		if(($depth !~ /all/) and ($depth<scalar(@total))){@select=(shuffle(@total))[0..$depth-1];}else{@select=@total;}

		open(NFSA,">./$depth.$permu/$n/$bar.fasta") or die "Cannot open new $file\n";
		foreach my $seqname(@select){
			$P{$bar}{$seqname}{"seq"} = $B{$bar}{$seqname}{"seq"};
			$P{$bar}{$seqname}{"info"} = $B{$bar}{$seqname}{"info"};
			$P{$bar}{$seqname}{"length"} = length($B{$bar}{$seqname}{"seq"});
			print NFSA ">",$seqname,"\t",$P{$bar}{$seqname}{"info"},"\n",$P{$bar}{$seqname}{"seq"},"\n";
		}
		close NFSA;
	}

	my %G;
	my @bararray;
	foreach my $bar (keys %P){
		push(@bararray,$bar);
		foreach my $genus (keys %G){
			push(@{$G{$genus}},0);
		}
		foreach my $seqname (keys %{$P{$bar}}){
			my $bn = scalar(@bararray);
			my $genus = $P{$bar}{$seqname}{"info"};
	if($genus eq ""){print $bar,"\t1",$seqname,"1\n";}
			unless(exists $G{$genus}){
				my @new = (0) x $bn;
				$G{$genus} =\ @new;
			}
			$G{$genus}->[$bn-1] += 1;
		}
	}
	print "Genus Hash.....done!!\n";

	############    UNIFRAC OUTPUT    #############
	open(EXCEL,">./$depth.$permu/summary.$n.csv") or die "Failed to open Result File\n";
	print EXCEL "barcode";
	foreach my $bar(@bararray){print EXCEL ",",$bar;}
	foreach my $genus (sort keys %G){
		print EXCEL "\n$genus";
		foreach (@{$G{$genus}}){print EXCEL ",",$_;}
	}
	close EXCEL;
	print "Output Summary.....done!!\n";
}
system("cp all.1/summary.1.csv summary.genera.csv");

