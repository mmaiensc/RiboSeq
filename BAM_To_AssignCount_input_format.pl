#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;

my $input; 
my $output; 
my $mismatch;

GetOptions ("input=s" => \$input,
				"output=s" => \$output,
				"mismatch=i" => \$mismatch);

my $usage = "Usage: $0\n--input <SAM or BAM in> (required)\n--output <tab out> (required)\n--mismatch <Number of mismatch(s)>, <3. Default 0(required)\n\n";

if (!defined $input || !defined $output || !defined $mismatch ) {
	print $usage;
	exit();
}

if (defined $mismatch > 2) {
	print STDERR "Number of mismatches($mismatch) cannot be greater than 2 ($mismatch)\n";
	exit();
}

if( substr($input,length($input)-4,4) eq ".bam" ){
	# check for samtools
	my $samtools = `which samtools 2> /dev/null | wc -l`;
	chomp($samtools);
	if( $samtools != 1 ){die "Error in BAM_To_AssignCount_input_format.pl: samtools not found for input bam file";}
	open (INFILE,"samtools view $input|");
}
else{
	open (INFILE, "<$input");
}
open (OUTFILE, ">$output");

while(my $line=<INFILE>){
	chomp($line);
	if($line =~ /NM:i:(\d+)/&&$1<=$mismatch){
		if($line =~ /XS:i:*/){
			}
			else {
				my @datum = split("\t", $line);
				## length of the sequence.
				my $lenght_seq = length($datum[9]);
				if ($datum[1] == 16){
					## correct for negative strand.
					my $correct_cordinates = $datum[3] + $lenght_seq - 1;
					print OUTFILE ("$datum[9]\t1\t$lenght_seq\t$correct_cordinates\t-\t$mismatch\n");
				}elsif ($datum[1] == 0){
					print OUTFILE ("$datum[9]\t1\t$lenght_seq\t$datum[3]\t+\t$mismatch\n");
				}else{
					print "ERROR:$line\n" ;
				}
			}	
	}

}

close(INFILE);
close(OUTFILE);

