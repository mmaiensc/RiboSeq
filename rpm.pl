#!/usr/bin/perl

use warnings;
use strict;

my $geneCount_fname = $ARGV[0]; 
my $gene_Count_col =  $ARGV[1]; 
my $assigncount_F_fname = $ARGV[2];
my $assigncount_R_fname = $ARGV[3];
my $assigncount_col=$ARGV[4];
my $wig_track_txt=$ARGV[5];
my $OUT1=$ARGV[6];
my $OUT1_Wig=$ARGV[7];
my $OUT2=$ARGV[8];
my $OUT2_Wig=$ARGV[9];


$gene_Count_col = $gene_Count_col  - 1 ;
$assigncount_col= $assigncount_col -1 ;


open (FH, "< $geneCount_fname") or die("Unable to open file : $!");

### get summation for given column based on genecount data. 

my $total = 0;
while (my $line=<FH>){
	chomp($line);
	my @datum = split('\t', $line);
	$total = $total + $datum[$gene_Count_col];
}
close(FH);

rpm($assigncount_F_fname,$OUT1,$OUT1_Wig);
rpm($assigncount_R_fname,$OUT2,$OUT2_Wig);

sub rpm{
	my ($fname, $fout, $fwig) = @_;
	print("$fname\t$fout\t********************$fwig\n");
	open (FH1, "< $fname") or die("Unable to open file : $!");
	open (FOUT1, "> $fout") or die("Unable to open file : $!");
	open (FO2, "> $fwig")or die("Unable to open file : $!");


	print FO2 "track type=wiggle_0\n";
	print FO2 "variableStep chrom=";
	print FO2 "$wig_track_txt\n" ;
	
	while (my $line=<FH1>){
		chomp($line);
		my @datum = split('\t', $line);
		my $rpm =  sprintf("%.8f" , (($datum[$assigncount_col]/$total) *1000000) ) ;
		printf FOUT1 "$datum[0]\t$rpm\t$datum[1]\t$datum[2]\n";
		printf FO2 "$datum[0]\t$rpm\n";
		}
	close(FH1);
	close(FOUT1);
	close(FO2);
}


# open (FH1, "< $assigncount_F_fname");
# open (FH2, "< $assigncount_R_fname");
# 
# 
# open (FOUT1, "> $OUT1") ;
# 
# 
# 
# 
# open (FOUT2, "> $OUT2") ;
# 
# while (defined(my $val1 = <FH1>) and defined(my $val2 = <FH2>)) {
#     my @datum1 = split('\t', $val1);
#     my @datum2 = split('\t', $val2);
# 	my $rpm_F =  sprintf("%.8f" , (($datum1[$assigncount_col]/$total) *1000000) ) ;
# 	my $rpm_R =  sprintf("%.8f" , (($datum2[$assigncount_col]/$total) *1000000) );
# 	
# 	printf FOUT1 ("$datum1[0]\t$rpm_F\t$datum1[1]\t$datum1[2]");
# 	printf FOUT2 ("$datum2[0]\t$rpm_R\t$datum2[1]\t$datum2[2]");
# 	
# }
# 
# close(FH1);
# close(FH2);
# close(FOUT1);
# close(FOUT2);
