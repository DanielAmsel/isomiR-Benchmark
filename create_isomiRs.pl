#! /usr/bin/perl
use strict;
use warnings;
# Author : Daniel.Amsel@ime.fraunhofer.de
# Purpose: Create artificial microRNA isomiRS for testing
# Output : 7 multi-fasta files with microRNA variants


### Supported Types ###
# Templated 5' Additions 		--> using the precursor sequence : adding 1nt|2nt|3nt to 5' end of miRNA
# Templated 3' Additions		--> using the precursor sequence : addint 1nt|2nt|3nt to 3' end of miRNA
# 5' Deletions				--> deleting 1nt|2nt|3nt of 5' end of miRNA
# 3' Deletions				--> deleting 1nt|2nt|3nt of 3' end of miRNA
# 3' non-template Additions		--> adding A|AA|AAA ; T|TT|TTT ; G|GG|GGG ; C|CC|CCC to the 3' end
# SNP in seed region (1-8)		--> changing nt to the three other nts	(original : ATGC | mod : TTGC ; CTGC ; GTGC | mod : AAGC ; AGGC ; ACGC | ...)
# SNP in rest region (9-end)		--> same as seed


#### use miR-file
# mutate every sequence and add them to hash{seq} =\@(id1,..)
# gather ids by identical sequence

my $mirs_file		= $ARGV[0];		# file of mature miRs
my $precursor_file	= $ARGV[1];		# file of precursor miRs



if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help'|| $ARGV[0] eq '' || $ARGV[0] eq '--h' || $ARGV[0] eq '--help')
{
help();
exit;
}


################################################################################################
my %precursor_hash	= %{&parse_fasta($precursor_file)};
################################################################################################
# five prime add
my $five_prime_file		= "$mirs_file-FIVE_ADD.fa";
my $fp_c			= 0;
open(FIVE,">",$five_prime_file);
my %five_prime_hash		= %{&fp_add($mirs_file,\%precursor_hash)};
foreach(keys %five_prime_hash){
	$fp_c			+=1;
	my $five_prime_seq	= $_;
	my @five_prime_array	= @{$five_prime_hash{$five_prime_seq}};
	my %five_prime_array_hash;
	foreach(@five_prime_array){
		$five_prime_array_hash{$_}="";
	}
	foreach(keys %five_prime_array_hash){
		print FIVE ">$_|";
	}
	print FIVE "FA-$fp_c\n";
	print FIVE "$five_prime_seq\n";
}
close(FIVE);
################################################################################################
# three prime add
my $three_prime_file		= "$mirs_file-THREE_ADD.fa";
my $tp_c			= 0;
open(THREE,">",$three_prime_file);
my %three_prime_hash		= %{&tp_add($mirs_file,\%precursor_hash)};
foreach(keys %three_prime_hash){
	$tp_c			+=1;
	my $three_prime_seq	= $_;
	my @three_prime_array	= @{$three_prime_hash{$three_prime_seq}};
	my %three_prime_array_hash;
	foreach(@three_prime_array){
		$three_prime_array_hash{$_}="";
	}
	foreach(keys %three_prime_array_hash){
		print THREE ">$_|";
	}
	print THREE "TA-$tp_c\n";
	print THREE "$three_prime_seq\n";
}
close(THREE);
################################################################################################
# three prime deletion
my $three_sub_file		= "$mirs_file-THREE_DEL.fa";
my $ts_c			= 0;
open(THREE,">",$three_sub_file);
my %three_sub_hash		= %{&three_sub($mirs_file)};
foreach(keys %three_sub_hash){
	$ts_c			+=1;
	my $three_sub_seq	= $_;
	my @three_sub_array	= @{$three_sub_hash{$three_sub_seq}};
	my %three_sub_array_hash;
	foreach(@three_sub_array){
		$three_sub_array_hash{$_}="";
	}
	foreach(keys %three_sub_array_hash){
		print THREE ">$_|";
	}
	print THREE "TS-$ts_c\n";
	print THREE "$three_sub_seq\n";
}
close(THREE);

################################################################################################
# five prime deletion
my $five_sub_file		= "$mirs_file-FIVE_DEL.fa";
my $fs_c			= 0;
open(FIVE,">",$five_sub_file);
my %five_sub_hash		= %{&five_sub($mirs_file)};
foreach(keys %five_sub_hash){
	$fs_c			+=1;
	my $five_sub_seq	= $_;
	my @five_sub_array	= @{$five_sub_hash{$five_sub_seq}};
	my %five_sub_array_hash;
	foreach(@five_sub_array){
		$five_sub_array_hash{$_}="";
	}
	foreach(keys %five_sub_array_hash){
		print FIVE ">$_|";
	}
	print FIVE "FS-$fs_c\n";
	print FIVE "$five_sub_seq\n";
}
close(FIVE);

################################################################################################
# SNP Seed nt-wise mutation
my $snp_seed_file	= "$mirs_file-SNP-SEED.fa";
my $ss_c		= 0;
open(SNP,">",$snp_seed_file);
my %snp_seed_hash	= %{&snp_seed($mirs_file)};
foreach(keys %snp_seed_hash){
	$ss_c			+=1;
	my $snp_seed_seq	= $_;
	my @snp_seed_array	= @{$snp_seed_hash{$snp_seed_seq}};
	my %snp_seed_array_hash;
	foreach(@snp_seed_array){
		$snp_seed_array_hash{$_}="";
	}
	foreach(keys %snp_seed_array_hash){
		print SNP ">$_|";
	}
	print SNP "SS-$ss_c\n";
	print SNP "$snp_seed_seq\n";
}
close(SNP);
################################################################################################
# SNP pos 9 to end nt-wise mutation 
my $snp_rest_file	= "$mirs_file-SNP-REST.fa";
my $sr_c		= 0;
open(SNP,">",$snp_rest_file);
my %snp_rest_hash	= %{&snp_rest($mirs_file)};
foreach(keys %snp_rest_hash){
	$sr_c			+=1;
	my $snp_rest_seq	= $_;
	my @snp_rest_array	= @{$snp_rest_hash{$snp_rest_seq}};
	my %snp_rest_array_hash;
	foreach(@snp_rest_array){
		$snp_rest_array_hash{$_}="";
	}
	foreach(keys %snp_rest_array_hash){
		print SNP ">$_|";
	}
	print SNP "SR-$sr_c\n";
	print SNP "$snp_rest_seq\n";
}	
close(SNP);
################################################################################################
# non template nt addition
my $non_template_file		= "$mirs_file-NON-TEMPLATE.fa";
my $nt_c			= 0;
open(NTF,">",$non_template_file);
my %nt_file_hash		= %{&non_template($mirs_file)};
foreach(keys %nt_file_hash){
	$nt_c			+=1;
	my $nt_file_seq		= $_;
	my @nt_file_array	= @{$nt_file_hash{$nt_file_seq}};
	my %nt_file_array_hash;
	foreach(@nt_file_array){
		$nt_file_array_hash{$_}="";
	}
	foreach(keys %nt_file_array_hash){
		print NTF ">$_|";
	}
	print NTF "NT-$nt_c\n";
	print NTF "$nt_file_seq\n";
}
close(NTF);

################################################################################################
#######					 SUBROUTINES 					########
################################################################################################
# build a hash for precursor {tca-mir-1}=seq
sub parse_fasta{
	my $pf_file	= $_[0];
	my %pf_hash;
	my $pf_id;
	open(PF,"<",$pf_file) or die " Can't open $pf_file:\n$!\n";
	while(<PF>){
		chomp;
		my $pf_line	= $_;
		if(/^>/){
			my @pf_split	= split(" ",$pf_line);
			$pf_id		= $pf_split[0];
			$pf_id		=~s/^>//;
		}
		else{
			my $pf_seq	= $pf_line;
			$pf_seq		=~s/U/T/g;
			$pf_hash{$pf_id}=$pf_seq;
		}
	}
	close(PF);
	return(\%pf_hash);
}

# create snp seed
sub snp_seed{
	my $ss_file	= $_[0];
	my @ss_nt	= ('A','T','C','G');
	my %ss_hash;	#{seq}	=\@(id1,id2,..)
	my $ss_id;	# tca-miR-92c-5p
	open(SS,"<",$ss_file) or die " Can't open $pf_file:\n$!\n";
	while(<SS>){
		chomp;
		my $ss_line	= $_;
		if(/^>/){	# >tca-miR-92c-5p MIMAT0018786
			my @ss_id_split	= split(" ",$ss_line);
			$ss_id		= lc($ss_id_split[0]);
			$ss_id		=~s/>//; 
		}
		else{
			my $ss_seq	= $ss_line;
			$ss_seq		=~s/U/T/g;
			# mutate index 0-7
			# loop through seed sequence
			for(my $i=0; $i <= 7; $i++){	
				my $ss_mut_seq	= $ss_seq;
				my $ss_char	= substr($ss_seq,$i,1);
				# loop through all 4 bases to exchange 
				foreach(@ss_nt){
					next if ($ss_char eq $_);
					substr($ss_mut_seq,$i,1)=$_;					

					if(not exists $ss_hash{$ss_mut_seq}){
						my @ss_tmp1		= ($ss_id);
						$ss_hash{$ss_mut_seq}	= \@ss_tmp1;
					}
					else{
						my @ss_tmp2		= @{$ss_hash{$ss_mut_seq}};
						push(@ss_tmp2,$ss_id);
						$ss_hash{$ss_mut_seq}	= \@ss_tmp2;
					}
				}	
			}
		}
	}
	close(SS);	
	return(\%ss_hash);
}
# create snp rest
sub snp_rest{
	my $sr_file	= $_[0];
	my @sr_nt	= ('A','T','C','G');
	my %sr_hash;	#{seq} =\@(id1,id2,..)
	my $sr_id;
	open(SR,"<",$sr_file) or die " Can't open $pf_file:\n$!\n";
	while(<SR>){
		chomp;
		my $sr_line	= $_;
		if(/^>/){
			my @sr_id_split	= split(" ",$sr_line);
			$sr_id		= lc($sr_id_split[0]);
			$sr_id		=~s/>//;
		}
		else{
			my $sr_seq	= $sr_line;
			$sr_seq		=~s/U/T/g;
			my $sr_len	= length($sr_seq);
			for(my $j=8; $j<$sr_len;$j++){
				my $sr_mut_seq	= $sr_seq;
				my $sr_char	= substr($sr_seq,$j,1);
				foreach(@sr_nt){
					next if ($sr_char eq $_);
					substr($sr_mut_seq,$j,1)=$_;

					if(not exists $sr_hash{$sr_mut_seq}){
						my @sr_tmp1		= ($sr_id);
						$sr_hash{$sr_mut_seq}	= \@sr_tmp1;
					}
					else{
						my @sr_tmp2		= @{$sr_hash{$sr_mut_seq}};
						push(@sr_tmp2,$sr_id);
						$sr_hash{$sr_mut_seq}	= \@sr_tmp2;
					}
				}
			}
		}
	}
	close(SR);
	return(\%sr_hash);
}
# create 3' non-template-mod
sub non_template{
	my $nt_file	= $_[0];
	my @nt_nt	= ('A','T','C','G');
	my %nt_hash;
	my $nt_id;
	open(NT,"<",$nt_file) or die " Can't open $pf_file:\n$!\n";
	while(<NT>){
		chomp;
		my $nt_line	= $_;
		if(/^>/){
			my @nt_id_split	= split(" ",$nt_line);
			$nt_id		= lc($nt_id_split[0]);
			$nt_id		=~s/>//;
		}
		else{
			my $nt_seq	= $nt_line;
			$nt_seq		=~s/U/T/g;
			foreach(@nt_nt){
				my $nt_mut_seq	= $nt_seq;
				for ( my $l=0; $l < 3; $l++){
					$nt_mut_seq	.= $_;
					if(not exists $nt_hash{$nt_mut_seq}){
						my @nt_array1		= ($nt_id);
						$nt_hash{$nt_mut_seq}	= \@nt_array1;
					}				
					else{
						my @nt_array2		= @{$nt_hash{$nt_mut_seq}};
						push(@nt_array2,$nt_id);
						$nt_hash{$nt_mut_seq}	= \@nt_array2;
					}
				}
			}			
		}
	}
	close(NT);
	return(\%nt_hash);
}






##############################################################################################################################
# create 3' template-mod

# >tca-mir-3851a CGGACGAATTATGAATTTGGCGCCGAGTTGATTTGGAGTAGTACTATGGTTTTGGGTACTCACATAGTACATTCCGAATCGACTCGTTAACAAATTTAAATTCAACAA
# >tca-mir-3851a-1/2/3-t                  AGTTGATTTGGAGTAGTACTATGgtt
sub tp_add{
	my $tpa_file		= $_[0];
	my %tpa_prec_hash	= %{$_[1]};
	my %tpa_hash;
	my $tpa_id;
	my $tpa_idp;
	open(TPA,"<",$tpa_file) or die " Can't open $pf_file:\n$!\n";
	while(<TPA>){
		chomp;
		my $tpa_line	= $_;
		if(/^>/){
			my @tpa_split	= split(" ",$tpa_line);
			$tpa_id		= lc($tpa_split[0]);
			$tpa_id		=~s/>//;
			$tpa_idp	= $tpa_id;
			$tpa_id		=~s/-3p//;
			$tpa_id		=~s/-5p//;
		}
		else{
			my $tpa_seq	= $tpa_line;
			$tpa_seq	=~s/U/T/g;
			if(exists $tpa_prec_hash{$tpa_id}){
				my $tpa_preseq	= $tpa_prec_hash{$tpa_id};
				my $tpa_start	= index($tpa_preseq,$tpa_seq);
				my $tpa_len	= length($tpa_seq);
				my $tpa_pre_len	= length($tpa_preseq);
				for(my $o=1;$o<=3;$o++){
					my $tpa_addpos	= $tpa_start+$tpa_len;
					next if ($tpa_addpos+$o>$tpa_pre_len);
					my $tpa_addseq	= substr($tpa_preseq,$tpa_addpos,$o);
					my $tpa_newseq	= $tpa_seq.$tpa_addseq;
					if(not exists $tpa_hash{$tpa_newseq}){
						my @tpa_array1		= ($tpa_idp);
						$tpa_hash{$tpa_newseq}	= \@tpa_array1;
					}
					else{
						my @tpa_array2		= @{$tpa_hash{$tpa_newseq}};
						push(@tpa_array2,$tpa_idp);
						$tpa_hash{$tpa_newseq}	= \@tpa_array2;
					}
				}
			}
			else{
				foreach(keys %tpa_prec_hash){
					if($_ =~ /^$tpa_id/){
						my $tpa_preID_alt	= $_;
						my $tpa_preseq_alt	= $tpa_prec_hash{$tpa_preID_alt};
						my $tpa_start_alt	= index($tpa_preseq_alt,$tpa_seq);
						my $tpa_len_alt		= length($tpa_seq);
						my $tpa_pre_len_alt	= length($tpa_preseq_alt);
						for(my $p=1; $p <= 3 ; $p++){
							my $tpa_addpos_alt	= $tpa_start_alt+$tpa_len_alt;
							next if ($tpa_addpos_alt+$p>$tpa_pre_len_alt);
							my $tpa_addseq_alt	= substr($tpa_preseq_alt,$tpa_addpos_alt,$p);
							my $tpa_newseq_alt	= $tpa_seq.$tpa_addseq_alt;
							if(not exists $tpa_hash{$tpa_newseq_alt}){
								my @tpa_array3			= ($tpa_idp);
								$tpa_hash{$tpa_newseq_alt}	= \@tpa_array3;
							}
							else{
								my @tpa_array4			= @{$tpa_hash{$tpa_newseq_alt}};
								push(@tpa_array4,$tpa_idp);
								$tpa_hash{$tpa_newseq_alt}	= \@tpa_array4;
							}
						}
					}
				}
			}
		}	
	}
	close(TPA);	
	return(\%tpa_hash);
}
##############################################################################################################################
# create 5' template-mod

# >tca-mir-3851a CGGACGAATTATGAATTTGGCGCCGAGTTGATTTGGAGTAGTACTATGGTTTTGGGTACTCACATAGTACATTCCGAATCGACTCGTTAACAAATTTAAATTCAACAA
# >tca-mir-3851a-1/2/3-f               ccgAGTTGATTTGGAGTAGTACTATG
sub fp_add{
	my $fpa_file		= $_[0];	# the mature mir file
	my %fpa_prec_hash	= %{$_[1]};	# the precursor hash {id}=seq
#	my %fpa_mat_hash;			# mature mir hash{id}=seq
	my %fpa_hash;				# return_hash{$seq}=\@(id1,id2,..)
	my $fpa_id;
	my $fpa_idp;
	open(FPA,"<",$fpa_file) or die " Can't open $pf_file:\n$!\n";
	while(<FPA>){
		chomp;
		my $fpa_line	= $_;
		if(/^>/){			# trim ID 
			my @fpa_split	= split(" ",$fpa_line);
			$fpa_id		= lc($fpa_split[0]);
			$fpa_id		=~s/>//;
			$fpa_idp	= $fpa_id;
			$fpa_id		=~s/-3p//;
			$fpa_id		=~s/-5p//;
		}
		else{
			my $fpa_seq	= $fpa_line;
			$fpa_seq	=~s/U/T/g;
			# search seq in precursor (compare IDs)
			if(exists $fpa_prec_hash{$fpa_id}){
				my $fpa_preseq	= $fpa_prec_hash{$fpa_id};
				my $fpa_start	= index($fpa_preseq,$fpa_seq);	# start position of the fpa_seq in fpa_preseq
				# get nt one,two,three pos before and add it to the seq
				for(my $m = 1; $m <=3; $m++){
					my $fpa_addpos	= $fpa_start-$m;
					my $fpa_addseq	= substr($fpa_preseq,$fpa_addpos,$m);	# nucleotide(s) add
					my $fpa_newseq	= $fpa_addseq.$fpa_seq;
					if(not exists $fpa_hash{$fpa_newseq}){
						my @fpa_array1		= ($fpa_idp);
						$fpa_hash{$fpa_newseq}	= \@fpa_array1;
					}
					else{
						my @fpa_array2		= @{$fpa_hash{$fpa_newseq}};
						push(@fpa_array2,$fpa_idp);
						$fpa_hash{$fpa_newseq}	= \@fpa_array2;
					}
				}
			}
			else{ # some mature miRs can have two precursors (-1 or -2) 
				foreach(keys %fpa_prec_hash){
					if($_ =~ /^$fpa_id/){
						my $fpa_preID_alt	= $_;
						my $fpa_preseq_alt	= $fpa_prec_hash{$fpa_preID_alt};
						my $fpa_start_alt	= index($fpa_preseq_alt,$fpa_seq);
						for(my $n=1; $n <= 3; $n++){
							my $fpa_addpos_alt	= $fpa_start_alt-$n;
							my $fpa_addseq_alt	= substr($fpa_preseq_alt,$fpa_addpos_alt,$n);
							my $fpa_newseq_alt	= $fpa_addseq_alt.$fpa_seq;
							if(not exists $fpa_hash{$fpa_newseq_alt}){
								my @fpa_array3			= ($fpa_idp);
								$fpa_hash{$fpa_newseq_alt}	= \@fpa_array3;
							}
							else{
								my @fpa_array4			= @{$fpa_hash{$fpa_newseq_alt}};
								push(@fpa_array4,$fpa_idp);
								$fpa_hash{$fpa_newseq_alt}	= \@fpa_array4;
							}
						}
					}
				}
			}
		}
	}
	close(FPA);
	return(\%fpa_hash);
}

# create 5' shortened
sub five_sub{
        my $fs_file     = $_[0];
        my %fs_hash;
        my $fs_id;
        open(FS,"<",$fs_file) or die " Can't open $pf_file:\n$!\n";
        while(<FS>){
                chomp;
                my $fs_line     = $_;
                if(/^>/){
                        my @fs_id_split = split(" ",$fs_line);
                        $fs_id          = lc($fs_id_split[0]);
                        $fs_id          =~s/>//;
                }
                else{
                        my $fs_seq      = $fs_line;
                        $fs_seq         =~s/U/T/g;
			my $fs_mut_seq	= $fs_seq;
                        for ( my $r=0; $r < 3; $r++){
                        	$fs_mut_seq     =reverse($fs_mut_seq); # shorten
				chop($fs_mut_seq);
				$fs_mut_seq	=reverse($fs_mut_seq);
				
                                if(not exists $fs_hash{$fs_mut_seq}){
                                	my @fs_array1           = ($fs_id);
                                        $fs_hash{$fs_mut_seq}   = \@fs_array1;
                                }
                                else{
                                	my @fs_array2           = @{$fs_hash{$fs_mut_seq}};
                                        push(@fs_array2,$fs_id);
                                        $fs_hash{$fs_mut_seq}   = \@fs_array2;
                                }
                        }
                }
        }
        close(FS);
        return(\%fs_hash);
}

#create 3' shortened
sub three_sub{
        my $ts_file     = $_[0];
        my %ts_hash;
        my $ts_id;
        open(TS,"<",$ts_file) or die " Can't open $pf_file:\n$!\n";
        while(<TS>){
                chomp;
                my $ts_line     = $_;
                if(/^>/){
                        my @ts_id_split = split(" ",$ts_line);
                        $ts_id          = lc($ts_id_split[0]);
                        $ts_id          =~s/>//;
                }
                else{
                        my $ts_seq      = $ts_line;
                        $ts_seq         =~s/U/T/g;
                        my $ts_mut_seq  = $ts_seq;
                        for ( my $s=0; $s < 3; $s++){
                                chop($ts_mut_seq);
                                if(not exists $ts_hash{$ts_mut_seq}){
                                	my @ts_array1           = ($ts_id);
                                        $ts_hash{$ts_mut_seq}   = \@ts_array1;
                                }
                                else{
                               		my @ts_array2           = @{$ts_hash{$ts_mut_seq}};
                                        push(@ts_array2,$ts_id);
                                        $ts_hash{$ts_mut_seq}   = \@ts_array2;
                                }
                        }
                }
        }
        close(TS);
        return(\%ts_hash);
} 




sub help{
        print "###################################################################\n";
        print "./create_isomiRs.pl <mature_miRs.fasta> <precursor_mirs.fasta>\n";
        print "###################################################################\n";
        print "Contact: daniel.amsel\@ime.fraunhofer.de\n";
}
