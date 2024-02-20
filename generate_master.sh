#!/bin/bash
# run per-sample processing on an input fastq file for ribosomal processing
# input is a fastq file, and list of barcodes
# output are wig files, genes data, and pos/neg RPM data

set -e

# read in arguments
in1=""
in2=""
manifest1=""
manifest2=""
tar1=""
tar2=""
ref_genome=""
ref_annotation=""
out=""
outmaster=""

# options
# gene list cutoff
cutoff=100
# cutoffs for master file
cutoff1=5
cutoff2=5
# how to normalize codons
norm='segment'

clean="y"
path=""

# usage statement
usage="USAGE:
	REQUIRED:
	-i1 input sample name 1
	-i2 input sample name 2
	-m1 input manifest 1
	-m2 input manifest 2
	    NOTE: sample name 1 is the control, sample name 2 is the treatment
	-G  reference genome file (fasta)
	-A  reference annotation file
	    format: [ID] [start] [end] [strand] [name]
	-o  output file (tab-delimited text)

	OPTIONAL:
	-c  cutoff threshold for common_genes (default $cutoff)
	-c1 cutoff threshold for master file, sample 1 (default $cutoff1)
	-c2 cutoff threshold for master file, sample 2 (default $cutoff2)
	-norm do codon RPM relative to the whole gene, or segmented by first
	      10 codons, last 10 codons, and in between ('all' or 'segment', default $norm)

	-clean remove intermediate files (y or n, default $clean)

	-path  path to executable scripts, if not in PATH
"

# parse command line
for (( i=1; i<= $#; i++)); do
	j=$((i+1))
	arg="${!i}"
	val="${!j}"
	case $arg in
		"-i1")
			in1="$val"
			i=$j
			;;
		"-i2")
			in2="$val"
			i=$j
			;;
		"-m1")
			manifest1="$val"
			i=$j
			;;
		"-m2")
			manifest2="$val"
			i=$j
			;;
		"-G")
			ref_genome="$val"
			i=$j
			;;
		"-A")
			ref_annotation="$val"
			i=$j
			;;
		"-o")
			out="$val"
			i=$j
			;;
		"-c")
			cutoff="$val"
			i=$j
			;;
		"-c1")
			cutoff1="$val"
			i=$j
			;;
		"-c2")
			cutoff2="$val"
			i=$j
			;;
		"-norm")
			norm="$val"
			i=$j
			;;
		"-clean")
			clean="$val"
			i=$j
			;;
		"-path")
			path="$val"
			i=$j
			;;
	*)
		echo -e "Error: flag $arg not recognized.\n$usage"
		exit 1
	esac
done

# check required arguments
if [ "$in1" = "" ]; then
	echo -e "Error: set -i1\n$usage"
	exit 1
fi
if [ "$in2" = "" ]; then
	echo -e "Error: set -i2\n$usage"
	exit 1
fi
if [ "$manifest1" = "" ]; then
	echo -e "Error: set -m1\n$usage"
	exit 1
fi
if [ "$manifest2" = "" ]; then
	echo -e "Error: set -m2\n$usage"
	exit 1
fi
if [ "$out" = "" ]; then
	echo -e "Error: set -o\n$usage"
	exit 1
fi
if [ "$ref_genome" = "" ]; then
	echo -e "Error: set -G\n$usage"
	exit 1
fi
if [ "$ref_annotation" = "" ]; then
	echo -e "Error: set -A\n$usage"
	exit 1
fi
if [ "$norm" != "segment" ] && [ "$norm" != "all" ]; then
	echo -e "Error: set -norm to be 'segment' or 'all'\n$usage"
	exit 1
fi
if [ "$norm" = "segment" ]; then
	do_segment="y"
else
	do_segment="n"
fi

# check that manifests exist
if [ ! -f "$manifest1" ]; then
	echo "Error: manifest $manifest1 not found"
	exit 1
fi
if [ ! -f "$manifest2" ]; then
	echo "Error: manifest $manifest2 not found"
	exit 1
fi
# get names of tar1 and tar2 from manifests
tar1=$(head -1 $manifest1 | awk '{if($3=="tar.gz") print $1}')
tar2=$(head -1 $manifest2 | awk '{if($3=="tar.gz") print $1}')
if [ "$tar1" = "" ]; then
	echo -e "Error: can't find a tarball in $manifest1\n"
	exit 1
fi
if [ "$tar2" = "" ]; then
	echo -e "Error: can't find a tarball in $manifest2\n"
	exit 1
fi
# get the paths to the manifests
man1_dir=$(cd $(dirname $manifest1) && pwd -P)
man2_dir=$(cd $(dirname $manifest2) && pwd -P)
# if tarballs don't exist, add the manifest directory
if [ ! -f "$tar1" ]; then
	tar1="$man1_dir/$tar1"
fi
if [ ! -f "$tar2" ]; then
	tar2="$man2_dir/$tar2"
fi
# check again if we can find them
if [ ! -f "$tar1" ]; then
	echo -e "Error: can't find tarball $tar1";
	exit 1
fi
if [ ! -f "$tar2" ]; then
	echo -e "Error: can't find tarball $tar2";
	exit 1
fi

# check that the sample names exist in the manifests
check1=$(awk -F "\t" -v n="$in1" '{if($2==n) print $2}' < $manifest1 | uniq | wc -l)
if (( $check1 != 1 )); then
	echo "Error: sample name $in1 not found in $manifest1"
	exit 1
fi
check2=$(awk -F "\t" -v n="$in2" '{if($2==n) print $2}' < $manifest2 | uniq | wc -l)
if (( $check2 != 1 )); then
	echo "Error: sample name $in2 not found in $manifest2"
	exit 1
fi



#
# utility functions
#
# find executable for $1: either itself, or the $path/$1
function find {
	check=$(which $1 2> /dev/null | wc -l | awk '{print $1}')
	if [ "$check" = 1 ]; then
		echo "$1"
	elif [ "$path" != "" ]; then
		check=$(ls $path/$1 2> /dev/null | wc -l | awk '{print $1}')
		if [ "$check" = 1 ]; then
			echo "$path/$1"
		else
			(>&2 echo "Error: script $1 can't be found")
			exit 1
		fi
	else
		(>&2 echo "Error: script $1 can't be found")
		exit 1
	fi
}

# in manifest $1 and tarball $2, look for the file matching sample name $3 and type $4
# depends on output directory $tmp existing
# extract that file from the tarball, and return the name
function get_file {
	manifest=$1
	tarball=$2
	sample=$3
	type=$4
	# find the file name; only return the first match; check if there's more than 1??
	name=$(awk -F "\t" -v n="$sample" -v t="$type" '{if($2==n && $3==t) print $1}' < $manifest | head -1)
	# check that this file is actually in the tarball
	check=$(tar -tvf $tarball | grep -w $name | wc -l)
	if (( $check != 1 )); then
		(>&2 echo "Error: can't find file $name in tarball $tar, though it is in the manifest $manifest")
		exit 1
	fi
	# extract it
	tar -xf $tarball -C $tmp $name
	# return the name
	name="$tmp/$name"
	echo "$name"
}

#
# start the pipeline
#
# load dependencies
common_genes=$(find common_genes.pl)
codon_rpm=$(find codon_rpm.pl)
fold_change_master_flexible=$(find fold_change_master_flexible.pl)
samtools=$(find samtools)
final=$(find final.pl)

# get reference fasta name
gname=$(head -1 $ref_genome | cut -f1 -d ' ' | sed 's/>//')

# make output tmp folder
tmp="$out.tmp"
if [ ! -d "$out.tmp" ]; then
	mkdir $tmp
fi


#
# start processing
#
# common genes
echo -e "# Getting common genes"
genes1=$(get_file $manifest1 $tar1 $in1 "Gene_counts")
genes2=$(get_file $manifest2 $tar2 $in2 "Gene_counts")
$common_genes $genes1 $genes2 $cutoff 5 $tmp/genes.list

# codon rpm
echo -e "# Computing codon RPMs"
codon_pos1=$(get_file $manifest1 $tar1 $in1 "Gene_posRPM")
codon_neg1=$(get_file $manifest1 $tar1 $in1 "Gene_negRPM")
$codon_rpm $codon_pos1 $codon_neg1 $ref_annotation $tmp/genes.list $tmp/codonPosRPM1.txt $tmp/codonNegRPM1.txt $do_segment

codon_pos2=$(get_file $manifest2 $tar2 $in2 "Gene_posRPM")
codon_neg2=$(get_file $manifest2 $tar2 $in2 "Gene_negRPM")
$codon_rpm $codon_pos2 $codon_neg2 $ref_annotation $tmp/genes.list $tmp/codonPosRPM2.txt $tmp/codonNegRPM2.txt $do_segment

# make master file
echo -e "# Making master tables"
$fold_change_master_flexible $tmp/codonPosRPM1.txt $tmp/codonPosRPM2.txt $ref_annotation $ref_genome $gname $cutoff1 $cutoff2 9 F $tmp/MasterPos.txt $samtools 2> /dev/null
$fold_change_master_flexible $tmp/codonNegRPM1.txt $tmp/codonNegRPM2.txt $ref_annotation $ref_genome $gname $cutoff1 $cutoff2 9 R $tmp/MasterNeg.txt $samtools 2> /dev/null

# merge runs
echo -e "# Combining master tables for positive and negative"
$final $tmp/MasterPos.txt $tmp/MasterNeg.txt $out

# compress intermediates
if [ "$clean" != "y" ]; then	
	tar -zcf $out.tmp.tar.gz $tmp
fi
rm -r $tmp

exit 0
