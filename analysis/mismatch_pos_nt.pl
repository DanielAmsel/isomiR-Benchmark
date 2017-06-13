#! /usr/bin/perl
use strict;
use warnings;

my $result_file	= $ARGV[0];
my $condition	= $ARGV[1];
my $threshold	= 0;

my %result_hash;	
$result_hash{'freq'}	= 0;
#$result_hash{'AT'}	= 0;	# A overwrites T
#$result_hash{'AC'}      = 0;
#$result_hash{'AG'}      = 0;
#$result_hash{'TA'}      = 0;
#$result_hash{'TG'}      = 0;
#$result_hash{'TC'}      = 0;
#$result_hash{'CG'}      = 0;
#$result_hash{'CA'}      = 0;
#$result_hash{'CT'}      = 0;
#$result_hash{'GA'}      = 0;
#$result_hash{'GC'}      = 0;
#$result_hash{'GT'}      = 0;




open(IN,"<",$result_file);
while(<IN>){
	chomp;
	my $result_line		= $_;
	my @result_split	= split("\t",$result_line);
	my $result_freq_total	= $result_split[2];
	my $result_mism		= $result_split[6];
	my $result_ambigue	= $result_split[14];

	next if ($.==1);
	next if ($result_freq_total < $threshold);	# minimal read abundancy
	my $result_freq		= $result_freq_total/$result_ambigue;
	$result_hash{'freq'}	+= $result_freq;
	next if ($result_mism eq "0");
	if(not exists $result_hash{$result_mism}){
		$result_hash{$result_mism} = $result_freq;
	}
	else{
		$result_hash{$result_mism}+=$result_freq;
	}
	
}
close(IN);

print "position;mm;cpm\n";
foreach(keys %result_hash){
	my $id	= $_;
	next if ($id eq 'freq');
	my $cpm	= $result_hash{$_}/($result_hash{'freq'}*1000000);
	my @split_id	= split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $id);
	my @mutation	= split("",$split_id[1]);
	print "$condition;$split_id[0];$mutation[1]>$mutation[0];$cpm\n";
}

