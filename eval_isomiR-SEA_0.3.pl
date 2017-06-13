#! /usr/bin/perl
use strict;
use warnings;
# isomiRID is sensitive for ...-1 ...-2 precursor distinguishing (tca-mir-3851a-1, tca-mir-3851a-2,...)
# adapted this by cutting -\d away from the end in the result file  
my $unique_file         = $ARGV[0];
my $ambigue_file        = $ARGV[1];
my $mir_file            = $ARGV[2];
my $species             = $ARGV[3];



my %len_hash;           # {15}=\@(TP,FP,FN);
my @len_array           = (15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30);
foreach(@len_array){
        my @la_tmp_array= (0,0,0);
        $len_hash{$_}   = \@la_tmp_array;
}



my %mir_hash		= %{&parse_fasta($mir_file)};


my %unique_hash		= %{&parse_result($unique_file)};
my %ambigue_hash	= %{&parse_result($ambigue_file)};

my $unique_keys		= scalar( keys %unique_hash);
my $ambigue_keys	= scalar( keys %ambigue_hash);
my %result_hash		= (%unique_hash, %ambigue_hash);
my $result_keys		= scalar(keys %result_hash);

# should not happen, if isomiR-SEA is cool..
if($ambigue_keys+$unique_keys != $result_keys){
	print STDERR "Uneven HASHES ! $ambigue_keys+$unique_keys != $result_keys\n";
	exit;
}

my $tp	= 0;	# true  positive
my $fp	= 0;	# false positive
my $fn	= 0;	# false negative


# check for false negatives
foreach(keys %mir_hash){
	my $mh_mir_seq	= $_;
	my @mh_mir_ids	= @{$mir_hash{$mh_mir_seq}};
	# length
        my $mh_mir_len          = length($mh_mir_seq);  # length of sequence 
        my @mh_mir_ary          = @{$len_hash{$mh_mir_len}};    # tp fp fn
        my $mh_mir_tp           = $mh_mir_ary[0];
        my $mh_mir_fp           = $mh_mir_ary[1];
        my $mh_mir_fn           = $mh_mir_ary[2];

	if(not exists $result_hash{$mh_mir_seq}){
		$fn	+= scalar(@mh_mir_ids);
                $mh_mir_fn      += scalar(@mh_mir_ids);	# length
	}
	else{
		my @mh_result_ids	= @{$result_hash{$mh_mir_seq}};
		my %mh_result_id_hash	= map{$_=>1} @mh_result_ids;
		foreach(@mh_mir_ids){
			if(not exists $mh_result_id_hash{$_}){
				$fn	+=1;
                                $mh_mir_fn      += 1;
			}
		}
	}
        @mh_mir_ary             = ($mh_mir_tp,$mh_mir_fp,$mh_mir_fn);
        $len_hash{$mh_mir_len}  = \@mh_mir_ary;
}
# check for true positives and false positives
foreach(keys %result_hash){
	my $rh_result_seq	= $_;
	my @rh_result_ids	= @{$result_hash{$rh_result_seq}};
	# length
        my $rh_mir_len          = length($rh_result_seq);
        my @rh_mir_ary          = @{$len_hash{$rh_mir_len}};    # tp fp fn
        my $rh_mir_tp           = $rh_mir_ary[0];
        my $rh_mir_fp           = $rh_mir_ary[1];
        my $rh_mir_fn           = $rh_mir_ary[2];

	my @rh_mir_ids		= @{$mir_hash{$rh_result_seq}};
	
	my %rh_mir_id_hash	= map{$_=>1} @rh_mir_ids;
	foreach(@rh_result_ids){
		my $rh_result_id	= $_;
		if(not exists $rh_mir_id_hash{$rh_result_id}){
			$fp		+= 1;
                        $rh_mir_fp      += 1;
		}
		else{
			$tp		+= 1;
                        $rh_mir_tp      += 1;
		}
	}
        @rh_mir_ary             = ($rh_mir_tp,$rh_mir_fp,$rh_mir_fn);
        $len_hash{$rh_mir_len}  = \@rh_mir_ary;
}




print "TP\tFP\tFN\n";
print "$tp\t$fp\t$fn\n";

print "LEN\tTP\tFP\tFN\n";
foreach(@len_array){
        my $lh_key      = $_;
        my @lh_ary      = @{$len_hash{$lh_key}};
        print "$lh_key\t$lh_ary[0]\t$lh_ary[1]\t$lh_ary[2]\n";
}


sub parse_result{
	my $pr_file	= $_[0];
	my %pr_hash;
	
	open(PR,"<",$pr_file);
	while(<PR>){
		chomp;
		next if ($.==1);
		my $pr_line	= $_;
		my @pr_split	= split('\|',$pr_line); # split only works for ambigue file
		my $pr_read_seq;
		foreach(@pr_split){
			my $pr_result		= $_;
			$pr_result		=~s/^\t//;
			my @pr_result_array	= split("\t",$pr_result); # tag_index       tag_sequence    tag_quality     #count_tags     mirna_id        mirna_name      mirna_seq       begin_ungapped_mirna ....
			if(scalar(@pr_result_array)>=14){
				shift(@pr_result_array);
				$pr_read_seq	= shift(@pr_result_array);
				$pr_read_seq		=~s/U/T/g;
				shift(@pr_result_array);
				shift(@pr_result_array);
			}
			my $pr_mir_header	= $pr_result_array[1];
			my @pr_mir_header_array	= split(" ",$pr_mir_header);	# >tca-miR-10-5p MIMAT0008357 Tribolium castaneum miR-10-5p
			my $pr_mir_id		= lc($pr_mir_header_array[0]);
			$pr_mir_id		=~s/>//g;
			if(not exists $pr_hash{$pr_read_seq}){
				my @pr_tmp1		= ($pr_mir_id);
				$pr_hash{$pr_read_seq} 	= \@pr_tmp1;
			}
			else{
				my @pr_tmp2		= @{$pr_hash{$pr_read_seq}};
				push(@pr_tmp2,$pr_mir_id);
				$pr_hash{$pr_read_seq} 	=\@pr_tmp2;
			}	
		}	
		
	}
	close(PR);
	return(\%pr_hash);
}

# trim >
# split at |
sub parse_fasta{
        my $pf_file     = $_[0];
        my %pf_hash;
        my $pf_id;
        open(PF,"<",$pf_file);
        while(<PF>){
                chomp;
                my $pf_line     = $_;
                if(/^>/){
                        $pf_id  = lc($pf_line);
                        $pf_id  =~s/>//g;
                }
                else{
                        my $pf_seq      = $pf_line;
                        my @pf_id_split = split('\|',$pf_id);
                        pop(@pf_id_split);	# remove tailing counter
                        if(not exists $pf_hash{$pf_seq}){
				# possible double IDs are filtered out
				my %pf_id_split_hash;
				foreach(@pf_id_split){
					$pf_id_split_hash{$_}="";
				}
				my @pf_id_split_unique	= keys %pf_id_split_hash;
                                $pf_hash{$pf_seq}       = \@pf_id_split_unique;
                        }
                        else{
                                print STDERR "Sequence already in HASH... $pf_seq\n";
                        }
                }
        }
        close(PF);
        return(\%pf_hash);
}

