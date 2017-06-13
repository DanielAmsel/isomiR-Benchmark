#! /usr/bin/perl
use strict;
use warnings;

my $input	= $ARGV[0];
my $condition	= $ARGV[1];


# cpm ; condition ; len(-3,-2,-1,1,2,3)
my %result_hash_3;
my %result_hash_5;
$result_hash_3{'freq'}		= 0;
$result_hash_3{-3}		= 0;
$result_hash_3{-2}              = 0;
$result_hash_3{-1}              = 0;
$result_hash_3{1}               = 0;
$result_hash_3{2}               = 0;
$result_hash_3{3}               = 0;

$result_hash_5{'freq'}          = 0;
$result_hash_5{-3}              = 0;
$result_hash_5{-2}              = 0;
$result_hash_5{-1}              = 0;
$result_hash_5{1}               = 0;
$result_hash_5{2}               = 0;
$result_hash_5{3}               = 0;



my %five_del_seq;	# {id} = count
my %five_add_seq;	# {id} = count
my %three_del_seq;	# {id} = count
my %three_add_seq;	# {id} = count
my %all_types_seq;	# {id} = count <- count all unique reads

open(IN,"<",$input);
while(<IN>){
        chomp;
        next if ($.==1);
        my $line                = $_;
        my @split               = split("\t",$line);
	my $seq			= $split[0];
	my $id			= $split[1];
        my $freq                = $split[2];
        my $t5            	= $split[8];	# t5
	my $t3			= $split[9];	# t3
	my $ambigue		= $split[14];
	my $normalized_freq	= $freq/$ambigue;
        $result_hash_3{'freq'}  += $normalized_freq;
	$result_hash_5{'freq'}	+= $normalized_freq;


	$all_types_seq{$id}	= "";

        if ($t5 ne 0){
		my $t5_len	= length($t5);	
		if($t5 eq uc $t5){	# add
			$result_hash_5{"$t5_len"}+=$normalized_freq;
			if (not exists $five_add_seq{$id}){
				my @fas_array1		= ();
				push(@fas_array1,$t5_len);
				$five_add_seq{$id}	= \@fas_array1;
			}
			else{
				my @fas_array2		= @{$five_add_seq{$id}};
				push(@fas_array2,$t5_len);
				$five_add_seq{$id}	= \@fas_array2;
			}
		}
		elsif($t5 eq lc $t5){	# del
			$result_hash_5{"-$t5_len"}+=$normalized_freq;
                        if (not exists $five_del_seq{$id}){
                                my @fds_array1          = ();
                                push(@fds_array1,"-$t5_len");
                                $five_del_seq{$id}      = \@fds_array1;
                        }
                        else{
                                my @fds_array2          = @{$five_del_seq{$id}};
                                push(@fds_array2,"-$t5_len");
                                $five_del_seq{$id}      = \@fds_array2;
                        }

		}
	}
	elsif($t3 ne 0){
		my $t3_len	= length($t3);
		if($t3 eq uc $t3){	# add
			$result_hash_3{"$t3_len"}+=$normalized_freq;


                        if (not exists $three_add_seq{$id}){
                                my @tas_array1          = ();
                                push(@tas_array1,$t3_len);
                                $three_add_seq{$id}      = \@tas_array1;
                        }
                        else{
                                my @tas_array2          = @{$three_add_seq{$id}};
                                push(@tas_array2,$t3_len);
                                $three_add_seq{$id}      = \@tas_array2;
                        }

#			$three_add_seq{$id} = $t3_len;
		}
		elsif($t3 eq lc $t3){	# del
			$result_hash_3{"-$t3_len"}+=$normalized_freq;



                        if (not exists $three_del_seq{$id}){
                                my @tds_array1          = ();
                                push(@tds_array1,"-$t3_len");
                                $three_del_seq{$id}      = \@tds_array1;
                        }
                        else{
                                my @tds_array2          = @{$three_del_seq{$id}};
                                push(@tds_array2,"-$t3_len");
                                $three_del_seq{$id}      = \@tds_array2;
                        }




#			$three_del_seq{$id} = "-$t3_len";
		}
	}
	
}
close(IN);

# %five_add_id{$value_5}=nr_count
my %five_id;
$five_id{-3}	= 0;
$five_id{-2}    = 0;
$five_id{-1}    = 0;
$five_id{1}     = 0;
$five_id{2} 	= 0;
$five_id{3}	= 0;

foreach(keys %five_add_seq){
	my $fas_key	= $_;
	my @fas_array	= @{$five_add_seq{$fas_key}};
	my $fas_size	= scalar(@fas_array);
	foreach(@fas_array){	# (-1, 3 , ..) # possibly more than one position for each read
		$five_id{$_} += (1/$fas_size);
	}
}

#my %five_del_id;
#$five_del_id{-3}        = 0;
#$five_del_id{-2}        = 0;
#$five_del_id{-1}        = 0;
#$five_del_id{1}         = 0;     
#$five_del_id{2}         = 0;
#$five_del_id{3}         = 0;

foreach(keys %five_del_seq){
        my $fds_key     = $_;
        my @fds_array   = @{$five_del_seq{$fds_key}};
        my $fds_size    = scalar(@fds_array);
        foreach(@fds_array){    # (-1, 3 , ..) # possibly more than one position for each read
                $five_id{$_} += (1/$fds_size);
        }
}


my %three_id;
$three_id{-3}        = 0;
$three_id{-2}        = 0;
$three_id{-1}        = 0;
$three_id{1}         = 0;     
$three_id{2}         = 0;
$three_id{3}         = 0;

foreach(keys %three_add_seq){
        my $tas_key     = $_;
        my @tas_array   = @{$three_add_seq{$tas_key}};
        my $tas_size    = scalar(@tas_array);
#	print "$tas_key\t@tas_array\n";
	foreach(@tas_array){
		$three_id{$_} += (1/$tas_size);
	}
}

#my %three_del_id;
#$three_del_id{-3}        = 0;
#$three_del_id{-2}        = 0;
#$three_del_id{-1}        = 0;
#$three_del_id{1}         = 0;
#$three_del_id{2}         = 0;
#$three_del_id{3}         = 0;

foreach(keys %three_del_seq){
        my $tds_key     = $_;
        my @tds_array   = @{$three_del_seq{$tds_key}};
        my $tds_size    = scalar(@tds_array);
        foreach(@tds_array){
                $three_id{$_} += (1/$tds_size);
        }
}





#value : position -3 -2 -1 0 1 2 3
#pos 5' or 3'
print "condition;value;cpm;position;cpm_fill\n";
foreach(keys %result_hash_5){
	my $value_5 	= $_;
	next if ($value_5 eq 'freq');
	my $cpm_5	= $result_hash_5{$value_5}/$result_hash_5{'freq'}*1000000;	
	my $cpm_5_nr	= $five_id{$value_5}/$result_hash_5{'freq'}*1000000;
#	print "$value_5 : $result_hash_5{$value_5}\n";
	print "$condition;$value_5;$cpm_5;5';$cpm_5_nr\n";
}

foreach(keys %result_hash_3){
	my $value_3	= $_;
	next if ($value_3 eq 'freq');
	my $cpm_3	= $result_hash_3{$value_3}/$result_hash_3{'freq'}*1000000;
	my $cpm_3_nr	= $three_id{$value_3}/$result_hash_3{'freq'}*1000000;
#	print "$value_3 : $result_hash_3{$value_3} || $three_id{$value_3} / $result_hash_3{'freq'} * 1000000\n";
	print "$condition;$value_3;$cpm_3;3';$cpm_3_nr\n";
}

#print "3' $result_hash_3{'freq'} || 5' $result_hash_5{'freq'}\n";






