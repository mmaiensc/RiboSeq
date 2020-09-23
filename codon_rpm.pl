#!/usr/bin/perl

use warnings;
use strict;

# if(scalar(@ARGV) != 2){
#     print ("Usage: raw_density_file genepaires\n") ;
#     exit(1);
# }

my $rpm_forward = $ARGV[0];  ## has four columns : 31	0.0618346	0.333333333333	1 [ pos,rpm,raw_density,count ]
my $rpm_reverse = $ARGV[1];  ## has four columns : 31	0.0618346	0.333333333333	1 [ pos,rpm,raw_density,count ]
my $gene_annotation = $ARGV[2]; 
my $gene_pair = $ARGV[3];  
my $OUT1=$ARGV[4];
my $OUT2=$ARGV[5];


##Step1 : Read in RPM file for the forward strand
##Step2 : Read in RPM file for the reverse strand 



my %raw_density_neg;
my %raw_density_pos;

my %range_genes;

read_rmp_files($rpm_forward, "+");
read_rmp_files($rpm_reverse, "-");
get_geneList($gene_pair) ;


# my $count = keys %range_genes;
# print "$count range_genes\n" ;
#$count = keys %raw_density_neg;
#print "$count raw_density_neg\n" ;
#$count = keys %raw_density_pos;
#print "$count raw_density_pos\n" ;


open (FOUT1, "> $OUT1") ;
open (FOUT2, "> $OUT2") ;
## Header pos rpm rd count
print FOUT1 ("Gene_id\tPosition\trpm\traw_density\tcount\tcodon_relative_rpm\tcodon_rpm\tcodon_count\tcodon_raw_density\n");
## Header pos rpm rd count
print FOUT2 ("Gene_id\tPosition\trpm\traw_density\tcount\tcodon_relative_rpm\tcodon_rpm\tcodon_count\tcodon_raw_density\n");

get_codon();
close (FOUT1);
close (FOUT2);



sub read_rmp_files {
	my ($fname, $strand) = @_;
	### data eg : 22411 0.022748 0.25 1
	open (RD, "< $fname") or die("Unable to open file : $!");
	if($strand eq "-"){
		while(my $l=<RD>){
		chomp($l);
		my($pid, $rpm, $rd, $count) = split(" ", $l);
		my $value = $rpm. "\t" . $rd . "\t" .$count;
		if (exists $raw_density_neg{$pid}){
			print("raw_density_neg{$pid} Error -- same co-oridnate should not exist twice:\t$l\n");
		}		
		else {
			$raw_density_neg{$pid} = $value ;
			#print("$pid\t\t\t\t$raw_density_neg{$pid}\n");
		}
		}
		close(RD);
	}elsif ($strand eq "+"){
		while(my $l=<RD>){
		chomp($l);
		my($pid, $rpm, $rd, $count) = split(" ", $l);
		my $value = $rpm. "\t" . $rd . "\t" .$count;
		if (exists $raw_density_pos{$pid}){
			print("raw_density_pos{$pid} Error -- same co-oridnate should not exist twice:\t$l\n");
		}		
		else {
			$raw_density_pos{$pid} = $value ;
			#print("$pid\t\t\t\t$raw_density_pos{$pid}\n");
		}
	}
		close(RD);
	
	}
}

sub get_geneList{

	## b3646;yicG	3816897	3817511	+	inner membrane protein, UPF0126 family
	## b3221;yhcH	3367039	3367500	-	conserved protein
	my $fname = $_[0];
	open (INFILE, "< $fname") or die("Unable to open file : $!");
	while(my $id=<INFILE>){
		chomp($id);
		$id =~ s/"//g;	# added in case someone decides to put quotes in their gene names, because, you know (MMC 3/12/19)
		my $gene_details=`grep -w \"$id\" $gene_annotation` ;
		chomp($gene_details);
		my @datum = split("\t", $gene_details);	
		my $gene_length=0;
		my ($r1_start, $r1_end) ;
		my ($r2_start, $r2_end) ;
		my ($r3_start, $r3_end) ;
		my ($r1_range, $r2_range, $r3_range) ;
		if($datum[3] eq '-'){
			if($datum[1] < $datum[2]){
			my $swap = $datum[1];
			$datum[1] = $datum[2];
			$datum[2] = $swap ;
			}
			### Check what is the length of the gene
			if( ($datum[1] - $datum[2]) <= 60){
				## for short gene : whole gene is our gene range.
				$range_genes{$datum[0]} = $datum[3] . "\t" . $datum[2]  . "_" . $datum[1] ;
			}
			else{
			
				## First 30 bases 
				$r1_start = $datum[1] ; ## gene_start
				$r1_end = $datum[1] - 29 ; 
				$r1_range =  $r1_start . "_" . $r1_end ;
	
				### Last 30 bases 
				$r3_start = $datum[2] + 29 ;
				$r3_end =  $datum[2] ; 
				$r3_range =   $r3_start . "_"  . $r3_end ;
	
				$r2_start = $r1_end - 1;
				$r2_end = $r3_start + 1 ;
				$r2_range =  $r2_start . "_" . $r2_end;
				### updates the global hash
				$range_genes{$datum[0]} = $datum[3] . "\t" . $r1_range . "\t" . $r2_range . "\t" . $r3_range  ;
				#print("$datum[0]\t$datum[3]\t$r1_range\t$r2_range\t$r3_range\t$gene_length\n");
			}		
	}elsif($datum[3] eq '+'){
			if( ($datum[2] - $datum[1]) <= 60){
				$range_genes{$datum[0]} = $datum[3] ."\t".  $datum[1] .  "_" . $datum[2];
			}
			else {
	
				## First 30 bases 
				$r1_start = $datum[1] ; ## gene_start
				$r1_end = $datum[1] + 29 ; 
				$r1_range =  $r1_start . "_" . $r1_end ;
	
				$r3_start = $datum[2] - 29 ;
				$r3_end =  $datum[2] ; 
				$r3_range =   $r3_start . "_"  . $r3_end ;
	
				$r2_start = $r1_end + 1;
				$r2_end = $r3_start - 1 ;
				$r2_range =  $r2_start . "_" . $r2_end;
	
				$range_genes{$datum[0]} = $datum[3] . "\t" . $r1_range . "\t" . $r2_range . "\t" . $r3_range  ;
				#print("$datum[0]\t$datum[3]\t$r1_range\t$r2_range\t$r3_range\t$gene_length\n");
			}
		}
	}
	
		close(INFILE)
}


sub get_codon {
	foreach my $gid (keys %range_genes) {
		### First calculate 3 rpm values
		my @range = split("\t",$range_genes{$gid});  ### Strand plus three values in this array
		#print ("$gid ==== $range_genes{$gid}\n");
		for(my $i=1; $i<scalar(@range); $i++){
			my $total_rpm = 0;
			my ($start, $end) = split("_" , $range[$i]) ;
			
			if ($range[0] eq '-'){	
				#print ("$start\t$end\n");
				for(my $j=$start; $j>=$end; $j--){
					if(exists $raw_density_neg{$j}){
						my @datum = split("\t" , $raw_density_neg{$j});
						$total_rpm = $total_rpm + $datum[0];
						#print ("$j\t$raw_density_neg{$j}\n");
					}
				}
			} elsif($range[0] eq '+'){	
					for(my $j=$start; $j<=$end; $j++){
						if(exists $raw_density_pos{$j}){
							my @datum = split("\t" , $raw_density_pos{$j});
							$total_rpm = $total_rpm + $datum[0];
						}
					}
			} #end of else
			if($range[0] eq "+"){
				print FOUT1 ("$gid\t$start\t$end\t$total_rpm\n");
				get_codon_details_pos ($gid,$start, $end, $total_rpm);  
			}elsif ($range[0] eq "-"){
				print FOUT2 ("$gid\t$start\t$end\t$total_rpm\n");
				get_codon_details_neg ($gid,$start, $end, $total_rpm); 
				}
			} ## end of range
	} ## end of hash

} # end sub routine






sub get_codon_details_neg {

	my ($gene_id, $start, $end, $total_rpm) = @_;
	#print("$gene_id, $start, $end, $total_rpm\n");
	if($total_rpm == 0){
		print FOUT2 ("$gene_id\t$start\tTOTAL_RPM_ZERO\n");
	}else {
	## 3066102_3066073
	for(my $i=$start; $i>=$end; $i=$i-3){
		#print("i $i\n");
		my $total_codon_rpm = 0;
		my $total_codon_relative_rpm = 0 ;
		my $total_codon_count = 0;
		my $total_raw_density_count=0;
		my @codon;
		my $counter=0;
		if(exists $raw_density_neg{$i}){
		    $codon[0]="$gene_id\t$i\t$raw_density_neg{$i}";
		   ### data eg : 22411 0.022748 0.25 1 [ pos :      rpm rd count ]
			my @datum = split("\t" , $raw_density_neg{$i});
			#[ pos rpm rd count ]
			my $rpm_by_total_rpm = $datum[0] / $total_rpm ;
			$total_codon_relative_rpm = $total_codon_relative_rpm + $rpm_by_total_rpm ;
			$total_codon_count = $total_codon_count + $datum[2];
			$total_raw_density_count = $total_raw_density_count + $datum[1];
			$total_codon_rpm = $total_codon_rpm +  $datum[0];
		}else {
			$codon[0]="$gene_id\t$i\tNA\tNA\tNA";
		}
		#2nd element in the codon (triplet)
		$counter = $i - 1;
		if(exists $raw_density_neg{$counter}){
			$codon[1]="$gene_id\t$counter\t$raw_density_neg{$counter}";
			my @datum = split("\t" , $raw_density_neg{$counter});
			my $rpm_by_total_rpm = $datum[0] / $total_rpm ;
			$total_codon_relative_rpm = $total_codon_relative_rpm + $rpm_by_total_rpm ;
			$total_codon_count = $total_codon_count + $datum[2];
			$total_raw_density_count = $total_raw_density_count + $datum[1];
			$total_codon_rpm = $total_codon_rpm +  $datum[0];
		}else {
			$codon[1]="$gene_id\t$counter\tNA\tNA\tNA";
		}
		#3rd element in the codon (triplet)
		$counter = $i - 2;
		if(exists $raw_density_neg{$counter}){
			$codon[2]="$gene_id\t$counter\t$raw_density_neg{$counter}";
			my @datum = split("\t" , $raw_density_neg{$counter});
			my $rpm_by_total_rpm = $datum[0] / $total_rpm ;
			$total_codon_relative_rpm = $total_codon_relative_rpm + $rpm_by_total_rpm ;
			$total_codon_count = $total_codon_count + $datum[2];
			$total_raw_density_count = $total_raw_density_count + $datum[1];
			$total_codon_rpm = $total_codon_rpm +  $datum[0];
		}else {
			$codon[2]="$gene_id\t$counter\tNA\tNA\tNA";
		}
		
		if( ($total_codon_relative_rpm==0) && ($total_codon_count==0) ){
		    	print FOUT2 "$codon[0]\n$codon[1]\n$codon[2]\n" ;
		    }else {
		    ##print("Gene_id\tPosition\trpm\traw_density\tcount\tcodon_relative_rpm\tcodon_rpm\tcodon_count\tcodon_raw_density\n");
				print FOUT2 "$codon[0]\t$total_codon_relative_rpm\t$total_codon_rpm\t$total_codon_count\t$total_raw_density_count\n$codon[1]\n$codon[2]\n" ;

			}	
	}
	
	}	
	
}




sub get_codon_details_pos {

	my ($gene_id, $start, $end, $total_rpm) = @_;
	#print("$gene_id, $start, $end, $total_rpm\n");
	if($total_rpm == 0){
		print FOUT1 ("$gene_id\t$start\tTOTAL_RPM_ZERO\n");
	}else {
	for(my $i=$start; $i<=$end; $i=$i+3){
		my $total_codon_rpm = 0;
		my $total_codon_relative_rpm = 0 ;
		my $total_codon_count = 0;
		my $total_raw_density_count=0;
		my @codon; ### codon array represents : one codon : i.e. 3 Nucleotides i.e. codon[0] === first nucleotide 
		my $counter=0;
		if(exists $raw_density_pos{$i}){
		    $codon[0]="$gene_id\t$i\t$raw_density_pos{$i}";
		   ### data eg : 22411 0.022748 0.25 1 [ pos :      rpm rd count ]
			my @datum = split("\t" , $raw_density_pos{$i});
			#[ pos rpm rd count ]
			my $rpm_by_total_rpm = $datum[0] / $total_rpm ;
			$total_codon_relative_rpm = $total_codon_relative_rpm + $rpm_by_total_rpm ;
			$total_codon_count = $total_codon_count + $datum[2];
			$total_raw_density_count = $total_raw_density_count + $datum[1];
			$total_codon_rpm = $total_codon_rpm +  $datum[0];
		}else {
			$codon[0]="$gene_id\t$i\tNA\tNA\tNA";
		}
		#2nd element in the codon (triplet)
		$counter = $i + 1;
		if(exists $raw_density_pos{$counter}){
			$codon[1]="$gene_id\t$counter\t$raw_density_pos{$counter}";
			my @datum = split("\t" , $raw_density_pos{$counter});
			my $rpm_by_total_rpm = $datum[0] / $total_rpm ;
			$total_codon_relative_rpm = $total_codon_relative_rpm + $rpm_by_total_rpm ;
			$total_codon_count = $total_codon_count + $datum[2];
			$total_raw_density_count = $total_raw_density_count + $datum[1];
			$total_codon_rpm = $total_codon_rpm +  $datum[0];
		}else {
			$codon[1]="$gene_id\t$counter\tNA\tNA\tNA";
		}
		#3rd element in the codon (triplet)
		$counter = $i + 2;
		if(exists $raw_density_pos{$counter}){
			$codon[2]="$gene_id\t$counter\t$raw_density_pos{$counter}";
			my @datum = split("\t" , $raw_density_pos{$counter});
			my $rpm_by_total_rpm = $datum[0] / $total_rpm ;
			$total_codon_relative_rpm = $total_codon_relative_rpm + $rpm_by_total_rpm ;
			$total_codon_count = $total_codon_count + $datum[2];
			$total_raw_density_count = $total_raw_density_count + $datum[1];
			$total_codon_rpm = $total_codon_rpm +  $datum[0];
		}else {
			$codon[2]="$gene_id\t$counter\tNA\tNA\tNA";
		}
		
		if( ($total_codon_relative_rpm==0) && ($total_codon_count==0) ){
		    	print FOUT1 "$codon[0]\n$codon[1]\n$codon[2]\n" ;	
		    }else {
		    ##print("Gene_id\tPosition\trpm\traw_density\tcount\tcodon_relative_rpm\tcodon_rpm\tcodon_count\tcodon_raw_density\n");
				#print "$codon[0]\n$codon[1]\n$codon[2]\n" ;
			print FOUT1 "$codon[0]\t$total_codon_relative_rpm\t$total_codon_rpm\t$total_codon_count\t$total_raw_density_count\n$codon[1]\n$codon[2]\n" ;
			}	
	}
	
	}
}








