#! /usr/bin/perl
use strict;
use warnings;

my $result_folder       = $ARGV[0];
my $mature_mir          = $ARGV[1]; # fasta file of mature mirs to have the total number
my $threshold           = 0;

opendir DIR, $result_folder;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR;

my %mature_hash;        # {mir}=0
open(MIR,"<",$mature_mir);
while(<MIR>){
        chomp;
        next unless (/^>/);
        my $mir_line    = $_;
        my @mir_split   = split(" ",$mir_line);
        my $mirID       = $mir_split[0];
        $mirID          =~s/>//;
        $mature_hash{$mirID}=0;
}
close(MIR);

my @sorted_mirs = sort { lc($a) cmp lc($b) } (keys %mature_hash);

print "-";
foreach(@sorted_mirs){
        print ";$_";
}
print "\n";


foreach(@files){
        my $file        = $_;
        my $path        = "$result_folder"."$file";
        my %result_hash;
        open(IN,"<",$path);
        while(<IN>){
                chomp;
                next if ($.==1);
                my $result_line         = $_;
                my @result_split        = split("\t",$result_line);
                next unless ($result_split[6] eq 0);
                next unless ($result_split[7] eq 0);
                next unless ($result_split[8] eq 0);
                next unless ($result_split[9] eq 0);
                next if ($result_split[2] < $threshold);

                my $result_mir          = $result_split[3]; # tca-miR-90-5p
                my $freq_total          = $result_split[2];
                my $ambigue             = $result_split[14];
                my $freq_cor            = $freq_total/$ambigue;

                $result_hash{'freq'}    += $freq_cor;
                if(not exists $result_hash{$result_mir}){
                        $result_hash{$result_mir} = $freq_cor;
                }
                else{
                        $result_hash{$result_mir}+= $freq_cor;
                }

        }
        close(IN);
        my $sort_freq   = $result_hash{'freq'};
        print "$file";
        foreach(@sorted_mirs){
                my $sort_mir    = $_;
                if (exists $result_hash{$sort_mir}){
                        my $result_freq = $result_hash{$sort_mir};
                        my $cpm         = $result_freq/($sort_freq*1000000);
                        print ";$cpm";
                }
        }
        print "\n";



}


