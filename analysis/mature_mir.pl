#! /usr/bin/perl
use strict;
use warnings;

# for heatmap plot
# e.g. last 4 conditions


my $result_file	= $ARGV[0];
my $condition	= $ARGV[1];
my $threshold	= 0;

my %result_hash;
$result_hash{'freq'}	= 0;

open(RESULT,"<",$result_file);
while(<RESULT>){
        chomp;
        my $result_line         = $_;
        next if ($.==1);
        my @result_split        = split("\t",$result_line);
        next unless ($result_split[6] eq 0);
        next unless ($result_split[7] eq 0);
        next unless ($result_split[8] eq 0);
        next unless ($result_split[9] eq 0);
	next if ($result_split[2] < $threshold);
	
	my $result_mir		= $result_split[3]; # tca-miR-90-5p
	my $freq_total		= $result_split[2];
	my $ambigue		= $result_split[14];
	my $freq_cor		= $freq_total/$ambigue;
	
	$result_hash{'freq'}	+= $freq_cor;
	if(not exists $result_hash{$result_mir}){
		$result_hash{$result_mir} = $freq_cor;
	}
	else{
		$result_hash{$result_mir}+= $freq_cor;
	}

}
close(RESULT);

foreach(keys %result_hash){
	my $key	= $_;
	next if ($key eq 'freq');
	my $cpm	= $result_hash{$key}/$result_hash{'freq'}*1000000;
	print "$cpm;$condition;$key\n";
}





