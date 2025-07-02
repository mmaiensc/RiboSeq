#!/usr/bin/perl
# fix end-of-line characters from an input stream

while( $line = <> ){
	chomp($line);
	@l2 = split('\r',$line);
	for($j=0; $j<= $#l2; $j++){
		print "$l2[$j]\n";
	}
}

