#!/usr/bin/perl

use warnings;
use strict;


my $f1= $ARGV[0];
my $f2= $ARGV[1];
my $OUT1=$ARGV[2];



unless (-e ($f1 ||  $f2) ) {
	print "Cannot find specified file" ;
}

open (FH1, "< $f1") or die("Unable to open file : $!");
open (FOUT1, "> $OUT1") ;

while(my $line=<FH1>){
	chomp($line);
	print FOUT1 "$line\n";
}

close(FH1);

my $skip=0;

open (FH1, "< $f2") or die("Unable to open file : $!");

while(my $line=<FH1>){
	chomp($line);
	if($skip==0){
		$skip=1;
		next;
	}
	print FOUT1 "$line\n";
}

close(FOUT1);
