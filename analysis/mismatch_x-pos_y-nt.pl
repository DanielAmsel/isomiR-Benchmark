#! /usr/bin/perl
use strict;
use warnings;

my $result_file	= $ARGV[0];
my $condition	= $ARGV[1];
my $threshold	= 0;

my %result_hash;

open(IN,"<",$result_file);
while(<IN>){
        chomp;
        my $result_line         = $_;
        my @result_split        = split("\t",$result_line);
        my $result_freq_total   = $result_split[2];
        my $result_mism         = $result_split[6];
        my $result_ambigue      = $result_split[14];

        next if ($.==1);
        next if ($result_freq_total < $threshold);      # minimal read abundancy
        my $result_freq         = $result_freq_total/$result_ambigue;
        $result_hash{'freq'}    += $result_freq;
        next if ($result_mism eq "0");
        if(not exists $result_hash{$result_mism}){
                $result_hash{$result_mism} = $result_freq;
        }
        else{
                $result_hash{$result_mism}+=$result_freq;
        }

}
close(IN);


open(A,">","MM_$condition-a.csv");
open(T,">","MM_$condition-t.csv");
open(G,">","MM_$condition-g.csv");
open(C,">","MM_$condition-c.csv");

my %muation_hash;	# {original}	= \@(modified,pos,cpm);
foreach(keys %result_hash){
        my $id  = $_;
        next if ($id eq 'freq');
	my $cpm = $result_hash{$_}/$result_hash{'freq'}*1000000;
	my @split_id    = split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $id);
        my @mutation    = split("",$split_id[1]);
	my $pos		= $split_id[0];
	my $original	= $mutation[1];
	my $modified	= $mutation[0];
	
	if($original eq "A"){
		print A "A;A>$modified;$pos;$cpm\n";
	}
	elsif($original eq "T"){
        	print T "T;T>$modified;$pos;$cpm\n";
        }
        elsif($original eq "G"){
                print G "G;G>$modified;$pos;$cpm\n";
        }
        elsif($original eq "C"){
                print C "C;C>$modified;$pos;$cpm\n";
        }


	

        
}


close(A);
close(T);
close(G);
close(C);
