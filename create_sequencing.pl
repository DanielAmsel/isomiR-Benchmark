#! /usr/bin/perl
use strict;
use warnings;
# Author : Daniel.Amsel@ime.fraunhofer.de
# Purpose: Create artificially sequenced microRNA reads
# Output : multi-fastq file 

# group multifasta file of microRNAs into length clusters.

# for each group, construct an artificial data set with 1000x coverage with ART
# if a length is not represented, no reads are created


# path to ART - adapt to your local system
my $art		= "/opt/art_bin_MountRainier/art_illumina";

my $mir_file	= $ARGV[0]; 	# fasta file of microRNAs
my $method	= $ARGV[1];	# HS20 MSv1
my @mir_len	= (17,18,19,20,21,22,23,24,25,26,27,28,29,30);



if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help'|| $ARGV[0] eq '' || $ARGV[0] eq '--h' || $ARGV[0] eq '--help')
{
help();
exit;
}




foreach(@mir_len){
	my $round_len		= $_;
	my $file_pattern	= "$round_len-miRs.fa";
	my $mir_id;
	open(OUT,">",$file_pattern);
	open(MIR,"<",$mir_file);
	while(<MIR>){
		chomp;
		my $mir_line	= $_;
		if(/^>/){
			$mir_id	= $_;
		}
		else{
			my $mir_seq	= $_;
			my $mir_len	= length($mir_seq);
			next unless ($mir_len == $round_len);
			print OUT "$mir_id\n$mir_seq\n";
		}
	}
	close(MIR);
	close(OUT);
	# run ART
	my $art_file		= "$round_len-$method-ART";
	my $run	= "$art -c 1000 -ss $method -i $file_pattern -l $round_len -o $art_file";
	system("$run");
}
# merge & clean 
system("cat *ART.fq > $method-all.fq");
system("cat *ART.aln > $method-all.aln");
system("rm -rf *ART.*");
system("rm -rf *-miRs.fa");


my $aln_file	= "$method-all.aln";
my $fq_file	= "$method-all.fq";
my $fa_file	= "$method-all.fa";

system("fastq_to_fasta -i $fq_file -o $fa_file");

my %aln_hash    = %{&aln_parser($aln_file)};
my %fa_hash     = %{&fasta_parser($fa_file)};



my %redundancy_hash;	# save simulated reads into hash --> unify reads

my $plus_reads		= "$method-all-plusStrand.fa";
open(PLUS,">",$plus_reads);
foreach(keys %fa_hash){
        my $fa_key      = $_;
        my $fa_seq      = $fa_hash{$fa_key};
        next unless (exists $aln_hash{$fa_key});
	print PLUS "$fa_key\n$fa_seq\n";		# all plus strand reads 
	if(not exists $redundancy_hash{$fa_seq}){
		$redundancy_hash{$fa_seq}=$fa_key;
	}
	else{
		next;
	}
}
close(PLUS);


my $unique_reads	= "$method-all-plusStrand-uniqueReads.fa";
open(OUT,">",$unique_reads);
foreach(keys %redundancy_hash){
	print OUT "$redundancy_hash{$_}\n$_\n";
}
close(OUT);





########################################################################################

sub fasta_parser{
        my $fp_file     = $_[0];
        my %fp_hash;
        my $fp_header;
        open(FP,"<",$fp_file);
        while(<FP>){
                chomp;
                my $fp_line     = $_;
                if($fp_line     =~/^>/){
                        my @fp_split            = split(" ",$fp_line);
                        $fp_header              = $fp_split[0];
                        $fp_hash{$fp_header}    = "";
                }
                else{
                        $fp_hash{$fp_header}    .= $fp_line;
                }
        }
        close(FP);
        return(\%fp_hash);
}



sub aln_parser{
        my $ap_file     = $_[0];
        my %ap_hash;

        open(AP,"<",$ap_file);
        while(<AP>){
                chomp;
                next unless (/^>/);
                my $ap_line     = $_;
                my @ap_split    = split("\t",$ap_line);
                next unless ($ap_split[3] eq "+");
                my $ap_key              = $ap_split[1];
                $ap_key                 = ">$ap_key";
                $ap_hash{$ap_key}="";
        }
        close(AP);
        return(\%ap_hash);
}

sub help{
	print "###################################################################\n";
	print "create_sequencing.pl <multi.fasta> <Method of ART>\n";
        print "###################################################################\n";
	print "Please note: You have to speicify the path of ART in the sourcecode\n"; 
	print "###################################################################\n";
	print "Contact: daniel.amsel\@ime.fraunhofer.de\n";
}







