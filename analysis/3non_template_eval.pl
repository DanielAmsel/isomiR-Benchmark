#! /usr/bin/perl
use strict;
use warnings;


# count all reads
# hash for all 3' non template additions
# {AAA} = 1234
# {AA}  = 23
# heatmap
# cpm condition nucleotide 

my %result_hash;

my $input	= $ARGV[0];
my $condition	= $ARGV[1];

open(IN,"<",$input);
while(<IN>){
	chomp;
	next if ($.==1);
	my $line		= $_;
	my @split		= split("\t",$line);
	my $freq		= $split[2];
	my $add			= $split[7];
	my $ambigue		= $split[14];
	my $cor_freq		= $freq/$ambigue;
	$result_hash{'freq'}	+= $cor_freq;	
	next unless ($add ne 0);
	if(not exists $result_hash{$add}){
		$result_hash{$add} 	= $cor_freq;
	}
	else{
		$result_hash{$add}     += $cor_freq;
	}
}
close(IN);


foreach(keys %result_hash){
	my $nt	= $_;
	my $cpm	= $result_hash{$nt} / $result_hash{'freq'}*1000000;
	print "$cpm;$condition;$nt\n";
}
