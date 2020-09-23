#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;
use Data::Dumper;

### given 2 files : using certain cutoff : get genes common to both files.


my $input1=$ARGV[0];
my $input2=$ARGV[1];
my $cutoff=$ARGV[2];
my $cutoff_column=$ARGV[3];
my $output=$ARGV[4];


# GetOptions ("input1=s" => \$input1,
# 				"input2=s" => \$input2,
# 				"cutoff=i" => \$cutoff,
# 				"cutoff_column=i" => \$cutoff_column,
# 				"output=s" => \$output);
# 
# my $usage = "Usage: $0\n--input1 <File 1 : IA tab in> (required)\n--input2 <File 2 : IA tab in> (required)\n--cutoff <minimum total count>. Recommended cutoff >=100(required)\n--cutoff_column <IA tab in> (required)\n--output <Output File Name> (required)\nWe assume that first column has gene names/IDs\n\n";

if (!defined $input1 || !defined $input2 || !defined $cutoff  || !defined $cutoff_column|| !defined $output ) {
	print "usage";
	exit();
}

if (defined $cutoff < 0) {
	print STDERR "Cutoff cannot be negative. Recommended cutoff 100\n";
	exit();
}



## based on the cutoff and cutoff column keep certain records.

my %f1_hash ;
my %f2_hash;
my $cc = $cutoff_column - 1 ;

open (FH1, "<$input1");

while(my $line=<FH1>){
	my @datum1 = split('\t', $line);
	if( $datum1[$cc] >= $cutoff ){
		$f1_hash{$datum1[0]} = 1;
	}
}
close(FH1);


open (FH2, "<$input2");

while(my $line=<FH2>){
	my @datum1 = split('\t', $line);
	if( $datum1[$cc] >= $cutoff ){
		$f2_hash{$datum1[0]} = 1;
	}
}

close(FH2);


# while (defined(my $val1 = <FH1>) and defined(my $val2 = <FH2>)) {
#     my @datum1 = split('\t', $val1);
#     my @datum2 = split('\t', $val2);
# 	if( $datum1[$cc] >= $cutoff ){
# 		$f1_hash{$datum1[0]} = 1;
# 	}
# 	if( $datum2[$cc] >= $cutoff ){
# 		$f2_hash{$datum1[0]} = 1;
# 	}
# }



open (FH2, "<$input2");
close(FH2);

# read($input1, \%f1_hash);
# read($input2, \%f2_hash);
# 
# sub read{
# 	my ($fname, $hash_ref) = @_;
# 	open (FH, "< $fname") or die("Unable to open file : $!");
# 	while (my $line=<FH>){
# 		chomp($line);
# 		my @datum = split('\t', $line);
# 		if( $datum[$cc] >= $cutoff ){
# 			$hash_ref->{$datum1[0]} = 1;
# 		}
# 	}
# 	close(FH);
# }


#my %intersected_hash = @f1_hash{grep exists $f2_hash{$_}, keys %f1_hash};
# foreach my $key (keys %f1_hash) {
#      print "$key\t$f1_hash{$key} \n";
# }
#     
#     
# foreach my $key (keys %f2_hash) {
#      print "$key\t$f2_hash{$key} \n";
# }

open (FO, "> $output") ;

for my $key (keys  %f1_hash) {
	if (exists $f2_hash{$key}) {
		print FO "$key\n" ;
	}
}

close(FO);
   
