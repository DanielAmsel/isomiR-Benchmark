#! /usr/bin/perl
use strict;
use warnings;


my $result_file		= $ARGV[0];
my $condition		= $ARGV[1];
my %result_hash;			# {mism}	= count

my $result_freq_ID;
my $result_mism_ID;
my $result_add_ID;
my $result_t5_ID;
my $result_t3_ID;

open(IN,"<",$result_file);	
while(<IN>){
	chomp;
	my $result_line		= $_;
	my @result_split	= split("\t",$result_line);
	my $result_freq_total	= $result_split[2];
	my $result_mism		= $result_split[6];
	my $result_add		= $result_split[7];
	my $result_t5		= $result_split[8];
	my $result_t3		= $result_split[9];
	my $result_ambigue	= $result_split[14];
	if ( $.== 1){
		$result_freq_ID			= $result_freq_total;
		$result_mism_ID			= $result_mism;
		$result_add_ID			= $result_add;
		$result_t5_ID			= $result_t5;
		$result_t3_ID			= $result_t3;

		$result_hash{$result_freq_ID}	= 0;
		$result_hash{$result_mism}	= 0;
		$result_hash{$result_add}	= 0;
		$result_hash{$result_t5}	= 0;
		$result_hash{$result_t3}	= 0;	
	}
	else{
		my $result_freq		= $result_freq_total / $result_ambigue;
		# subroutine
		$result_hash{$result_freq_ID}	+= $result_freq;
		$result_hash{$result_mism_ID}	+= &check_mutation($result_mism,$result_freq);
		$result_hash{$result_add_ID}	+= &check_mutation($result_add,$result_freq);
		$result_hash{$result_t5_ID}	+= &check_mutation($result_t5,$result_freq);
		$result_hash{$result_t3_ID}	+= &check_mutation($result_t3,$result_freq);
	}
}
close(IN);

foreach(keys %result_hash){
	next if($_ eq "freq");
	my $value	= $result_hash{$_}/$result_hash{"freq"}*1000000;
	print "$condition;$_;$value\n";
}





sub check_mutation{
	my $cm_text	= $_[0];	# 0 or something else, like AAA
	my $cm_freq	= $_[1];	# count of reads
	
	if($cm_text eq 0){
		$cm_freq	= 0;
	}
	return($cm_freq);
}
