#! usr/bin/perl
#! use strict; use warnings;

$file = $ARGV[0];

############    Establish Hash %B->barcode->seqname->length, genus, seq    #############
my %S;
#open(MTAX, "<./$file.final.unique.taxonomy") or die "Failed to open Mothur Taxonomy File\n";
open(RTAX, "<./$file.final.rdp.taxonomy.80") or die "Failed to open RDP output File\n";
open(FFSA, "<./$file.final.redundant.fasta") or die "Failed to open Final Fasta File\n";
open(OFSA, "<./$file.shhh.trim.fasta") or die "Failed to open shhh Fasta File\n";
open(GRP, "<./$file.final.groups") or die "Failed to open Groups File\n";
#open(NAM, "<./$file.final.names") or die "Failed to open Name File\n";
############    Read Group File    ############
local $/ = "\n";
while(<GRP>){
	chomp;
	my ($seqname,$bar) = split(/\t/);
	$S{$seqname}{"barcode"} = $bar;
}
close GRP;

############    Read Original Fasta File    ############
my %all;
local $/ = ">";
$_ = <OFSA>;
while(<OFSA>){
        chomp;
        my @array = split(/\n+/);
        my ($seqname,$junk) = split(/\s+/,shift(@array));
        my $seq = join("",@array);
#       $seq=~ s/\s|-|\.//g;
        $all{$seqname}{"oseq"} = $seq;
}
close OFSA;

############    Read RDP Taxonomy File    ############
local $/ = "\n";
while(<RTAX>){
	chomp;
	my ($seqname,$info) = split(/\t/);
#	$info =~ s/(\(|\);|;)/\t/g;
#       $info =~ s/"//g;
	local $/ = "\t";
	chomp($info);
#       my ($domain,$dscore) = split(/\(|\)/,$temp[0]);
#       my ($genus,$gscore) = split(/\(|\)/,$temp[-1]);
        $all{$seqname}{"oinfo"} = $info;
}
close RTAX;

############    Read Mothur Taxonomy File    ############
#local $/ = "\n";
#while(<MTAX>){
#        chomp;
#        my ($seqname,$info) = split(/\t/);
#        $info =~ s/(\(|\);|;)/\t/g;
#        local $/ = "\t";
#        chomp($info);
#        $S{$seqname}{"finfo"} = $info;
#        $S{$seqname}{"oinfo"} = $RDP{$seqname}{"oinfo"};
#}
#close MTAX;

############    Read Final Fasta File    ############
local $/ = ">";
$_ = <FFSA>;
while(<FFSA>){
	chomp;
	my @array = split(/\n+/);
	my ($seqname,$junk) = split(/\s+/,shift(@array));
	my $seq = join("",@array);
#	$seq=~ s/\s|-|\.//g;
	$S{$seqname}{"fseq"} = $seq;
	$S{$seqname}{"oseq"} = $all{$seqname}{"oseq"};
	$S{$seqname}{"oinfo"} = $all{$seqname}{"oinfo"};
	$S{$seqname}{"finfo"} = $all{$seqname}{"oinfo"};
}
close FFSA;

############    Read Names File    ############
#local $/ = "\n";
#while(<NAM>){
#	chomp;
#	my ($uniq,$string) = split(/\t/);
#	my @other = split(/,/,$string);
#	foreach my $seqname (@other){
#		$S{$seqname}{"finfo"} = $S{$uniq}{"finfo"};
#		$S{$seqname}{"oinfo"} = $RDP{$seqname}{"oinfo"};
#	}
#}
#close NAM;

print "Read in Files ..... done!!\n";

############    Establish Hash %B->barcode->seqname->length, genus, seq    #############
my %B;
foreach my $seqname (keys %S){
	my $bar = $S{$seqname}{"barcode"};

	$B{$bar}{$seqname}{"fseq"} = $S{$seqname}{"fseq"};
	$B{$bar}{$seqname}{"oseq"} = $S{$seqname}{"oseq"};
	$B{$bar}{$seqname}{"finfo"} = $S{$seqname}{"finfo"};
	$B{$bar}{$seqname}{"oinfo"} = $S{$seqname}{"oinfo"};
#	$B{$bar}{$seqname}{"length"} = length($oseq);
#	print $info,"\n";
}
print "B->barcode->seqname->length, genus, seq .... done!!\n";

############    Generate seperate fasta files    #############
system("rm -rf barcode.fin.fasta");
system("mkdir barcode.fin.fasta");
system("rm -rf barcode.shhh.fasta");
system("mkdir barcode.shhh.fasta");
foreach my $bar(keys %B){
	open(FBAR,">./barcode.fin.fasta/$bar.fasta") or die "cannot open final/$bar.fasta\n";
	open(OBAR,">./barcode.shhh.fasta/$bar.fasta") or die "cannot open shhh/$bar.fasta\n";
	foreach my $seqname(keys %{$B{$bar}}){
		print FBAR ">",$seqname,"\t",$B{$bar}{$seqname}{"finfo"},"\n",$B{$bar}{$seqname}{"fseq"},"\n";
		print OBAR ">",$seqname,"\t",$B{$bar}{$seqname}{"oinfo"},"\n",$B{$bar}{$seqname}{"oseq"},"\n";
	}
	close OBAR;
	close FBAR;
}

