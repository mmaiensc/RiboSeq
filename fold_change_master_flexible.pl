#!/usr/bin/perl
#
use warnings;
use strict;

## Control File
## Treatment File
## annotation file
## reference file
## cutoff for common
## chromosome_id

### possible issues : due to cut-off uneven number of genes 

if(scalar(@ARGV) != 11){
    print ("Usage: CHANGE ME ids_file file1 file2\n") ;
    exit(1);
}


my $codon_control = $ARGV[0];  
my $codon_treatment = $ARGV[1]; 
my $gene_annotation = $ARGV[2]; 
my $ref_file = $ARGV[3];
my $chr_name =   $ARGV[4];
my $cut_off_C = $ARGV[5];
my $cut_off_T = $ARGV[6];
my $cut_off_column = $ARGV[7];
my $strand = $ARGV[8];
my $OUT1=$ARGV[9];
my $samtools=$ARGV[$#ARGV];

my %control;
my %treatment;
my %ann_hash;


read_codon_details($codon_control, "C", $cut_off_C);
read_codon_details($codon_treatment, "T", $cut_off_T);
read_ann($gene_annotation);


open (FOUT1, "> $OUT1") ;
print FOUT1 ("GeneID_Treated\tPosition\tcodon_relative_rpm\tcodon_rpm\tcodon_count\tcodon_raw_density\t");
print FOUT1 ("GeneID_WT\tPosition[\]\tcodon_relative_rpm\tcodon_rpm\tcodon_count\tcodon_raw_density\t");
print FOUT1 ("fd_codon_relative_rpm\tfd_codon_rpm\tfd_codon_count\tfd_codon_raw_density\t");
print FOUT1 ("start_AA_cordinate\tend_AA_cordinate\tseq\trev_comp\tAA\tgene_start\tgene_end\n");

## rel_AA_start\trel_AA_end\tcase\n");
### Now get the intersection



for my $key (keys  %treatment) {
	if (exists $control{$key}) {
		my $fd=fold_change($treatment{$key}, $control{$key});
		#print("$record1\t$record2\t$fd\n");
		print FOUT1 ($fd) ;
	}
}
close (FOUT1);

sub fold_change{
	my ($r1, $r2) = @_ ;
	chomp($r1);
	chomp($r2);
	my @rec1 = split("\t", $r1); ## treated
	my @rec2 = split("\t", $r2); ## wt
	##Gene_id__Position	codon_relative_rpm	codon_rpm	codon_count	codon_raw_density
	##b2951;yggS__3093122	0.733858475642554	3.161472	149	35.4260073261
	##Gene_id	Position	rpm	raw_density	count	codon_relative_rpm	codon_rpm	codon_count	codon_raw_density
	
	my $fd_codon_relative_rpm = 0;
	my $fd_codon_rpm = 0;
	my $fd_codon_count	= 0;
	my $fd_codon_raw_density = 0;
	
	if($rec2[5]==0){
		$fd_codon_relative_rpm =  "NA" ;
	} else {
		$fd_codon_relative_rpm = $rec1[5] / $rec2[5] ;
		}
		
	if($rec2[6]==0){
		$fd_codon_rpm = "NA" ;
	} else {
		$fd_codon_rpm = $rec1[6] / $rec2[6] ;
	}
	
	if($rec2[7]==0){
		$fd_codon_count = "NA" ;
	} else {
		$fd_codon_count = $rec1[7] / $rec2[7] ;
	}
	
	if($rec2[8] == 0){
		$fd_codon_raw_density = "NA" ;
	}
	else {
		$fd_codon_raw_density = $rec1[8] / $rec2[8] ;
	}
	
	if( ($rec1[0] ne $rec2[0]) &&  ($rec1[1] != $rec2[1]) ){ #Geneid and position should be same for both of them
		print "Error : This should never happen :\t$rec1[0]\t$rec2[0]\n";
		exit;
	}
	my $rec1_all = $rec1[0] ."\t" . $rec1[1] . "\t" . $rec1[5] ."\t" . $rec1[6] . "\t" . $rec1[7] . "\t" . $rec1[8] ;
	my $rec2_all = $rec2[0] ."\t" . $rec2[1] . "\t" . $rec2[5] ."\t" . $rec2[6] . "\t" . $rec2[7] . "\t" . $rec2[8] ;
	my $fd_all = $fd_codon_relative_rpm . "\t" . $fd_codon_rpm . "\t" . $fd_codon_count . "\t" . $fd_codon_raw_density ;
	
	## get AA Seq
	my $AA = AAseq11($rec1[0], $rec1[1]) ;  ## pass gene name and position
	
	my $all = $rec1_all . "\t" . $rec2_all . "\t" . $fd_all . "\t" . $AA . "\n" ;
	#my $all = $r1 . "\t" . $r2 . "\t" . $fd_all ;
	
	return $all;
}


sub read_codon_details {
	my ($fname, $group, $cut_off) = @_;
	
	my $col_cut_off = $cut_off_column - 1;
	open (CD, "< $fname") or die("Unable to open file : $!");
	my $skip_header = 0;
	if($group eq "C"){
		while(my $l=<CD>){
		chomp($l);
		if($skip_header == 0){
			$skip_header++;
			next;
		}
		my @datum = split("\t", $l);
		if((scalar(@datum) >= 9)  && ($datum[$col_cut_off] >= $cut_off)){
			my $cid = $datum[0] . "__" . $datum[1] ;  ## geneName__position
			if (exists $control{$cid}){
				print("Control:$cid\t$control{$cid} Error -- same co-oridnate should not exist twice:\t$l\n");
			}		
			else {
				$control{$cid} = $l ;
				#print("$pid\t\t\t\t$control{$pid}\n");
			}
		}
		
		}
		close(RD);
		#print ("READ $fname\n");
	}elsif ($group eq "T"){
		while(my $l=<CD>){
		chomp($l);
		if($skip_header == 0){
			$skip_header++;
			next;
		}
		my @datum = split("\t", $l);
		if((scalar(@datum) >= 9)  && ($datum[$col_cut_off] >= $cut_off)){
			my $cid = $datum[0] . "__" . $datum[1] ;  ## geneName__position
			#print ("$cid\n");
			if (exists $treatment{$cid}){
				print("Treatment:$cid\t$treatment{$cid} Error -- same co-oridnate should not exist twice:\t$l\n");
			}		
			else {
				$treatment{$cid} = $l ;
				#print("$pid\t\t\t\t$control{$pid}\n");
			}
		}
	
	}
	close(RD);
	#print ("READ $fname\n");
	}
	
}




sub read_ann {
	my $annotation_fname = $_[0];
	open (FH, "< $annotation_fname") or die("Unable to open file : $!");
	while (my $line=<FH>){
	chomp($line);
	my @datum = split('\t', $line);
	##b3646;yicG	3816897	3817511	+	inner membrane protein, UPF0126 family
	if($datum[3] eq "-"){
		if($datum[1] < $datum[2]){
			my $swap = $datum[1];
			$datum[1] = $datum[2];
			$datum[2] = $swap ;
		}
	}
	my $value = $datum[1] . "\t" .$datum[2];
	$ann_hash{$datum[0]} = $value;
 }
}



sub AAseq11 {
	my ($gene_id, $location) = @_ ;
	#print ("$gene_id\t$location$strand\n");
	my ($gene_start,$gene_end) = split('\t' , $ann_hash{$gene_id});
	my $seq ;
	my $cordinates;
	my $rel_AA_start ;
	my $rel_AA_end ;
	my $case;
	
	##Krishna " While extracting the subsequence from the START codon, 
			##convert genes that START with GTG, TTG, CTG, ATT to ATG and 
			##then translate it to amino acid sequence."
	
	if($strand eq "F"){
		my $subtract = 27 ;
		my $add = 11;
		my $start_AA_cordinate = $location - $subtract ;
		my $end_AA_cordinate = $location + $add ;
		
		## if < start of gene,  < end of the gene : gene_start to Position+5
		## if < start of the gene , > end of the gene : gene_start to Gene_edn
		# if > gene start , > gene end  : pos -27 till end of the gene
		# default
		
		if(($start_AA_cordinate  < $gene_start ) && ( $end_AA_cordinate < $gene_end)) { 
			$seq = get_seq($gene_start, $end_AA_cordinate);
			## Check if seq begins with : 
			$seq =~ s/^(GTG|TTG|CTG|ATT)/ATG/ ;
			$cordinates = "$gene_start\t$end_AA_cordinate" ;
			$rel_AA_start = abs(  ($gene_start - $location) /3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 );
			$case = "C1_LT";
		}elsif(($start_AA_cordinate  < $gene_start ) && ( $end_AA_cordinate > $gene_end)){
			$seq = get_seq($gene_start, $gene_end);
			$cordinates = "$gene_start\t$gene_end" ;
			$rel_AA_start = abs(  ($gene_start - $location) /3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 );
			$case = "C2_GT_GE";
		}elsif(($start_AA_cordinate  > $gene_start ) && ( $end_AA_cordinate > $gene_end)){
			$seq = get_seq($start_AA_cordinate, $gene_end);
			$cordinates = "$start_AA_cordinate\t$gene_end" ;
			$rel_AA_start = abs(  ($gene_start - $location) /3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 ); 
			$case = "C3_GT";
		}else {
			$seq = get_seq($start_AA_cordinate, $end_AA_cordinate);
			$cordinates = "$start_AA_cordinate\t$end_AA_cordinate";
			$rel_AA_start = abs ( (($gene_start - $location ) )/3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 ); 
			#print("**** $cordinates\t$seq\n");
			$case = "C4_NOR";
		}
		my $rev_comp  = "";
		my $AA = codon($seq);
		if( !defined $AA ){
			$AA="Krishna_no11AA_pos";
		}
		# my $AA_details =  $cordinates . "\t" . $seq . "\t" .  $rev_comp . "\t" . $AA . "\t" . $gene_start . "\t" . $gene_end . "\t" . $rel_AA_start . "\t" . $rel_AA_end . "\t" . $case ;
		my $AA_details =  $cordinates . "\t" . $seq . "\t" .  $rev_comp . "\t" . $AA . "\t" . $gene_start . "\t" . $gene_end  ;

		#print $AA_details;
		return $AA_details;
	}elsif($strand eq "R"){
	
		###b2577;yfiE	2713342	2712461	-
		my $subtract = 11 ;
		my $add = 27;
		my $start_AA_cordinate = $location + $add ;
		my $end_AA_cordinate = $location - $subtract ;
		
		
		#print("S/E\t$start_AA_cordinate\t$end_AA_cordinate \n");
		
		#Check if position + 27 is > Start of gene. If no, this can be done normally. 
		#If yes, sequence we need is from codon position - 5 until the start of the gene.
		
		#Check if position - 5 is < Stop of the gene. If no, 
		#this can be done normally. If yes, sequence we need is from the end of the gene until position + 27.
		
		if(($start_AA_cordinate  > $gene_start ) && ( $end_AA_cordinate > $gene_end)) { 
			$seq = get_seq($end_AA_cordinate, $gene_start);
			## Check if seq begins with 
			$cordinates = "$end_AA_cordinate\t$gene_start" ;
			$rel_AA_start = abs(  ($gene_start - $location) /3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 ); 
			$case = "C1_GT_GS";
			#print("**** $cordinates\t$seq\n");
		}elsif( ($start_AA_cordinate  > $gene_start ) && ($end_AA_cordinate < $gene_end) ){
			$seq = get_seq($gene_start,$gene_end) ;
			$cordinates = "$gene_start\t$gene_end" ;
			$rel_AA_start = abs( (($gene_start - $location) )/3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 );
			$case = "C2_WG";
			#print("**** $cordinates\t$seq\n");
		}elsif( ($start_AA_cordinate < $gene_start ) && ($end_AA_cordinate < $gene_end) ){
			$seq = get_seq($start_AA_cordinate,$gene_end) ;
			$cordinates = "$start_AA_cordinate\t$gene_end" ;
			$rel_AA_start = abs( (($gene_start - $location) )/3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 ); #abs ( (($location- $gene_end) + 1)/3 );
			$case = "C3_LT_GS_GE";
		}else{
			$seq = get_seq($start_AA_cordinate, $end_AA_cordinate);
			$cordinates = "$start_AA_cordinate\t$end_AA_cordinate";
			$rel_AA_start = abs ( (($gene_start - $location ) )/3 );
			$rel_AA_end = abs ((($location- ($gene_end +1)))/3 ); 
			#print("**** $cordinates\t$seq\n");
			$case = "C4_NOR";
		}
		my $rev_comp = reverse_complement_IUPAC ($seq);
		$rev_comp =~ s/^(GTG|TTG|CTG|ATT)/ATG/ ;
		my $AA = codon($rev_comp);
		if( !defined $AA ){
			$AA="Krishna_no11AA_neg";
		}
		
		
		### Find where these AA are located in the gene
		#my $location = get_loc($ann_hash{$gene_id})

		my $start_pos = ($location); ### remains same for negative strand
		#my $t1 = $gene_start - $start_pos  ;
		#my $t2 =  $start_pos - $gene_end ;	
		#my $rel_codon_t1 = $t1 /3 ;
		#my $rel_codon_t2 = $t2 /3 ;
		
		my $AA_details =  $cordinates . "\t" . $seq . "\t" .  $rev_comp . "\t" . $AA . "\t" . $gene_start . "\t" . $gene_end  ;
		#my $AA_details =  $cordinates . "\t" . $seq . "\t" .  $rev_comp . "\t" . $AA . "\t" . $gene_start . "\t" . $gene_end . "\t" . $rel_AA_start . "\t" . $rel_AA_end . "\t" . $case ;
		#print $AA_details;
		return $AA_details;
	}

}


sub get_seq{
	my $start_pos = $_[0];
	my $end_pos = $_[1];
	if ($start_pos > $end_pos){ ## ? >=
		my $temp_swap = $end_pos;
		$end_pos = 	$start_pos;
		$start_pos = $temp_swap
	}
	my @seq = `$samtools faidx $ref_file $chr_name:${start_pos}-${end_pos}` ;
	chomp(@seq);
	my $combine="" ;
	for(my $i=1; $i<scalar(@seq); $i++){
		$combine = $combine . $seq[$i] ;
		#print ("$i\t$seq[$i]\n");
	}
	return($combine);
}

##http://code.izzid.com/2011/08/25/How-to-reverse-complement-a-DNA-sequence-in-perl.html

sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}



## https://docs.google.com/file/d/0BxxKmyjbL1eVZjZwN1gwXzR1NG8/edit?pli=1

sub codon{
	my $DNAseq = $_[0];
	my $codon ;
	my $protein ;
	
	for(my $i=0;$i<(length($DNAseq)-2);$i+=3){
	$codon=substr($DNAseq,$i,3);
	$protein .= codon2aa($codon);
	}
	return $protein ;

}


sub codon2aa{
	my($codon)=@_;
	$codon=uc $codon;
	my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

	if(exists $g{$codon}){return $g{$codon};}
	elsif($codon=~/GC./i){return 'A';}
	elsif($codon=~/GG./i){return 'G';}
	elsif($codon=~/CC./i){return 'P';}
	elsif($codon=~/AC./i){return 'T';}
	elsif($codon=~/GT./i){return 'V';}
	elsif($codon=~/CG./i){return 'R';}
	elsif($codon=~/TC./i){return 'S';}
	else{ return('x');
		print "Bad codon \"$codon\"!!\n";
		}
}



