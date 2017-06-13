#! /usr/bin/perl
use strict;
use warnings;

my $result_file	= $ARGV[0];
my $condition	= $ARGV[1];
my $freq	= 0;
my $count	= 0;
open(RESULT,"<",$result_file);
while(<RESULT>){
	chomp;
	my $result_line		= $_;
	next if ($.==1);
	my @result_split	= split("\t",$result_line);
	my $ambigue 		= $result_split[14];
	my $total_count		= $result_split[2];
	my $normalized_count	= $total_count / $ambigue;
	$count 			+=$normalized_count;
	next unless ($result_split[6] eq 0);
	next unless ($result_split[7] eq 0);
	next unless ($result_split[8] eq 0);
	next unless ($result_split[9] eq 0);
#	print "$result_line\n";	
	$freq 		+= $result_split[2];
}
close(RESULT);
my $cpm	= $freq/$count*1000000;


print "$condition;mature;$cpm\n";
