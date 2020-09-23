#!/usr/bin/perl
# get barcode from end of read and demultiplex based on it
#

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use lib dirname (__FILE__);

#
# read params
#
my $list = "";
my $out = "";
my $fastq = "";
my $fastq2 = "";
my $flex = 0;
my $min = 0;
&read_command_line();

#
# read in barcodes and sample names
#
(my $barcodes, my $blen) = &read_list( $list, $flex );


# 
# read and write fastq files
#
&write_fastq( $barcodes, $blen, $fastq, $fastq2, $flex, $min );


###########################
# read command line params
sub read_command_line{
	my $options = "Usage:
	REQUIRED:
	-i  Barcode list. Expects 2 columns, tab-delimited:
	      [barcode]  [sample name]
	    Outputs will be [sample name]_R1/2.fastq.gz
	    If [sample name] is omitted, the barcode will be used for the sample name.
	    This script looks for barcode matches at the end of the read in R1.
	-r1 Input R1

	OPTIONS:
	-o  Output prefix (prepended to output names, default: none)
	-r2 Input R2 (if available)
	-f  Flexible barcode search: look within the entire read for a match to the barcodes,
	    rather than strictly at the end of the read.
	-m  Minimum read length (mainly if using flexible barcode filtering, default $min)

	    WARNING: This will increase the search time: if you have N barcodes, it will take
	    about N times longer to run in most circumstances.

	    NOTE: If multiple matches are found, this script will report the one that's
	    closest to the end of the read.
";

	GetOptions('i=s' => \$list,
		   'o=s' => \$out,
		   'r1=s' => \$fastq,
		   'r2=s' => \$fastq2,
		   'f' => \$flex,
		   'm=i' => \$min
		) || die "$options\n";

	# check parameter settings
	if( $list eq "" ){
		die "Error: provide input -i\n$options\n";
	}
	if( $fastq eq "" ){
		die "Error: provide input -r1\n$options\n";
	}
}

# read in list of barcodes
sub read_list{
	my $b = {};
	my $blen = {};
	my $names = {};
	my $f = $_[1];
	open(IN,"$_[0]") || die "Error: can't open $_[0]\n";
	while(my $line = <IN>){
		chomp($line);
		my @parts = split('\t',$line);
	
		my $barcode = $parts[0];
		# get sample name: either provided in 2nd column, or use the barcode itself
		my $newname;
		if( $#parts >= 1 ){
			$newname = &fix_string( $parts[1] );
			if( $newname ne $parts[1] ){
				print "Warning: name $parts[1] changed to $newname to avoid special characters\n";
			}
		}
		else{
			$newname = $barcode;
		}

		# add name to list
		if( exists $names->{$newname} ){
			die "Error: after fixing the strings in your sample names, we have duplicate groups:\n$names->{$newname}\n$parts[1]\n\nPlease make your sample names more unique\n";
		}
		$names->{$newname} = $newname;

		if( exists $b->{$barcode} ){
			die "Error: barcode $barcode is repeated in $_[0]\n";
		}
		$b->{$barcode}->{name} = $newname;
		$b->{$barcode}->{count} = 0;

		my $len = length($barcode);
		if( $f ){
			$len = 1;
		}
		$blen->{$len} = 1;
	}
	close(IN);
	return ($b, $blen);
}

# fix string $_[0] to be good in use for filenames
sub fix_string{
	my $s = $_[0];
	$s=~ s/[ ?#*<>{}()@|\/\[\]\\"';:]/_/g;
	return $s;
}

# write a fastq file for each barcode from $_[0], using barcode lengths from $_[1], streaming fastq files $_[1] and $_[2] (if present)
# do search with regex if $_[4] is set; check read length against $_[5]
sub write_fastq{
	my $b = $_[0];
	my $blen = $_[1];
	my $r1 = $_[2];
	my $r2 = $_[3];
	my $f = $_[4];
	my $m = $_[5];
	my $nomatch = 0;

	# list of barcode lengths
	my @bl = keys %{$blen};

	# check we have gzip; use zcat instead
	#my $gzip = `which gzip 2> /dev/null | wc -l`;
	my $gzip = `which zcat 2> /dev/null | wc -l`;
	chomp($gzip);
	if( $gzip == 0 ){
		#die "Error: gzip not found in your PATH\n";
		die "Error: zcat not found in your PATH\n";
	}

	# open a stream for each group
	foreach my $barcode (keys %{$b}){
		open($b->{$barcode}->{r1},"| gzip -c > $out$b->{$barcode}->{name}_R1.fastq.gz");
		if( $r2 ne "" ){
			open($b->{$barcode}->{r2},"| gzip -c > $out$b->{$barcode}->{name}_R2.fastq.gz");
		}
	}
	
	# open r1, checking if gzipped or not
	#if( substr($r1,length($r1)-3,3) eq ".gz" ){
	#	open(R1,"gzip -dc $r1 |") || die "Error: can't open pipe to $r1\n";
	#}
	#else{
	#	open(R1,"$r1") || die "Error: can't open $r1\n";
	#}
	# open with zcat -f (force) regardless
	open(R1,"zcat -f $r1 |") || die "Error: can't open pipe to $r1\n";

	# open r2 if it exists, checking if gzipped or not
	if( $r2 ne "" ){
		#if( substr($r2,length($r2)-3,3) eq ".gz" ){
		#	open(R2,"gzip -dc $r2 |") || die "Error: can't open pipe to $r2\n";
		#}
		#else{
		#	open(R2,"$r2") || die "Error: can't open $r2\n";
		#}
		# open with zcat -f (force) regardless
		open(R2,"zcat -f $r2 |") || die "Error: can't open pipe to $r2\n";
	}

	# stream r1 (and r2, if present), write reads as we find them
	while(my $line1 = <R1>){
		# finish reading R1
		my $line2 = <R1>;
		my $line3 = <R1>;
		my $line4 = <R1>;
		chomp($line1);
		chomp($line2);
		chomp($line3);
		chomp($line4);
		my @parts = split(' ',$line1);
		my $r1id = $parts[0];
		$r1id =~ s/\@//;
		chomp($r1id);
		my @parts2 = split(':',$parts[$#parts]);
		# the sequence contains the barcode
		my $barcode = $line2;
		# read in R2 (if present)
		my $r2line1;
		my $r2line2;
		my $r2line3;
		my $r2line4;
		my $r2id;
		if( $r2 ne "" ){
			$r2line1 = <R2>;
			$r2line2 = <R2>;
			$r2line3 = <R2>;
			$r2line4 = <R2>;
			chomp($r2line1);
			chomp($r2line2);
			chomp($r2line3);
			chomp($r2line4);
			@parts = split(' ',$r2line1);
			$r2id = $parts[0];
			$r2id =~ s/\@//;
			chomp($r2id);
			# compare read IDs
			if( $r2id ne $r1id ){
				die "Error: fastq files out of order\n";
			}
		}

		# check if this barcode belongs to any of our lists
		# check substring matches to each barcode length
		for(my $i=0; $i<= $#bl; $i++){
			# barcode is from the end of the read
			if( length($barcode) > $bl[$i] ){
				my $newbarcode = "";
				my $pos = -1;	# start position of the barcode for flex search
				# for flexible search, find the closest to the end of the sequence
				if( $f ){
					$newbarcode = $barcode;
					my $bc = "";
					my @ret;
					foreach my $k (keys %{$b}){
						@ret = ();
						# find positions of all matches
						while( $newbarcode =~ /$k/g ){
							push(@ret, $-[0]);
						}
						# use the last one (closest to the end of the read)
						if( $#ret >= 0 && $ret[$#ret] > $m ){
							if( $ret[$#ret] > $pos ){
								$bc = $k;
								$pos = $ret[$#ret];
							}
						}
					}
					$newbarcode = $bc;
				}
				else{
					$newbarcode = substr($barcode,length($barcode)-$bl[$i],$bl[$i]);
				}
				if( exists $b->{$newbarcode} ){
					my $fout = $b->{$newbarcode}->{r1};
					$b->{$newbarcode}->{count}++;
					# trim the barcode off of the read
					if( $f ){
						# flexible search: everything before the barcode
						$line2 = substr($line2,0,$pos);
						$line4 = substr($line4,0,$pos);
					}
					else{
						# regular search: length of the barcode
						$line2 = substr($line2,0,length($line2)-$bl[$i]);
						$line4 = substr($line4,0,length($line4)-$bl[$i]);
					}
					print $fout "$line1\n$line2\n$line3\n$line4\n";
					if( $r2 ne "" ){
						$fout = $b->{$newbarcode}->{r2};
						print $fout "$r2line1\n$r2line2\n$r2line3\n$r2line4\n";
					}
				}
				else{
					$nomatch++;
				}
			}
		}
	}
	close(R1);
	if( $r2 ne "" ){
		close(R2);
	}

	# close output streams
	print "Sample Name\tNumber of reads\n";
	foreach my $barcode (keys %{$b}){
		close($b->{$barcode}->{r1});
		if( $r2 ne "" ){
			close($b->{$barcode}->{r2});
		}
		print "$b->{$barcode}->{name}\t$b->{$barcode}->{count}\n";
	}
	print "Undetermined\t$nomatch\n";
}







