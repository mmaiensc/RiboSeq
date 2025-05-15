#!/bin/bash
# run per-sample processing on an input fastq file for ribosomal processing
# input is a fastq file, and list of barcodes
# output are wig files, genes data, and pos/neg RPM data

set -e

# read in arguments
fastq=""
barcodes=""
ref_genome=""
ref_filter=""
ref_annotation=""
out=""
out_tar=""
out_manifest=""
tmp_tar=""
tmp_manifest=""
path=""
minreads=1000
flex="n"

# hard trim bases
trim5p=2
trim3p=5
# bowtie2 seeds
filter_seed=10
genome_seed=18
# alignment mismatches
mismatch=2
# psite params
pmethod="3p"
psite_offset=15
min_len=24
max_len=46

# usage statement
usage="USAGE:
	REQUIRED:
	-f  fastq data (can be zipped)
	-b  barcodes list (tab-delimited text)
	    Format: 1 barcode per line. Optionally include a second column with the output sample name.

	    Set this value to 'SKIP' (instead of a file name) to skip demultiplexing. Then it will use
	    'sample' as the same name. (Actually, providing anything that's not a valid file name will
	    have the same effect)
	-G  reference genome file (fasta)
	-A  reference annotation file
	    format: [ID] [start] [end] [strand] [name]
	-o  output prefix. Will be four output files:
	    [prefix].tar.gz: tarball of primary output results
	    [prefix].manifest.txt: manifest of contents in the tarball
	    [prefix].tmp.tar.gz: tarball of intermediate files
	    [prefix].tmp.manifest.txt: manifest of contents in tmp tarball

	    - or, specify each output in turn -

	-otar    [name]
	-oman    [name]
	-otmptar [name]
	-otmpman [name]

	    NOTE: the first line of each manifest will be the name of the output tarball. Keep these 
	    files together to keep this reference correct.

	OPTIONAL:
	-F     reference filter file (fasta)
	-min   minimum number of reads to analyze a sample (default $minreads)
	-flex  do demultiplexing flexibly: search for barcode in entire read (y or n, default $flex)
	       this is much more computationally intensive

	-hard5 bases to trim from 5' end of read (default $trim5p)
	-hard3 bases to trim from 3' end of read (default $trim3p)

	-Fseed bowtie2 seed length for filtering alignment (default $filter_seed)
	-Gseed bowtie2 seed length for genome alignment (default $genome_seed)

	-mis   max mismatches allowed in alignment (default $mismatch)

	-Pmethod P-site detection method: '3p', '5p', or 'center' (default $pmethod)
		 3p: P-site is Poff bases from the 3' end of the read
		 5p: P-site is Poff bases from the 5' end of the read
		 center: Trim Poff bases from both ends of the reads,
		 	 assign remaining P-site weight evenly over
			 remaining methods.
			 Reads shorter than 2*Poff+1 will be discarded.
	-Poff  P-site offset (default $psite_offset)
	-Pmin  minimum read length for P-site offset (default $min_len)
	-Pmax  maximum read length for P-site offset (default $max_len)

	-path  path to executable scripts, if not in PATH
"

# parse command line
for (( i=1; i<= $#; i++)); do
	j=$((i+1))
	arg="${!i}"
	val="${!j}"
	case $arg in
		"-f")
			fastq="$val"
			i=$j
			;;
		"-b")
			barcodes="$val"
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
		"-otar")
			out_tar="$val"
			i=$j
			;;
		"-oman")
			out_manifest="$val"
			i=$j
			;;
		"-otmptar")
			tmp_tar="$val"
			i=$j
			;;
		"-otmpman")
			tmp_manifest="$val"
			i=$j
			;;
		"-F")
			ref_filter="$val"
			i=$j
			;;
		"-flex")
			flex="$val"
			i=$j
			;;
		"-min")
			minreads="$val"
			i=$j
			;;
		"-hard5")
			trim5p="$val"
			i=$j
			;;
		"-hard3")
			trim3p="$val"
			i=$j
			;;
		"-Fseed")
			filter_seed="$val"
			i=$j
			;;
		"-Gseed")
			genome_seed="$val"
			i=$j
			;;
		"-mis")
			mismatch="$val"
			i=$j
			;;
		"-Pmethod")
			pmethod="$val"
			i=$j
			;;
		"-Poff")
			psite_offset="$val"
			i=$j
			;;
		"-Pmin")
			min_len="$val"
			i=$j
			;;
		"-Pmax")
			max_len="$val"
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
if [ "$fastq" = "" ]; then
	echo -e "Error: set -f\n$usage"
	exit 1
fi
if [ "$barcodes" = "" ]; then
	echo -e "Error: set -b\n$usage"
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
# check if output $out was given
if [ "$out" != "" ]; then
	out_manifest="$out.manifest.txt"
	tmp_manifest="$out.tmp.manifest.txt"
	out_tar="$out.results.tar.gz"
	tmp_tar="$out.tmp.tar.gz"
fi
# check that output names are defined
if [ "$out_tar" = "" ]; then
	echo -e "Error: set -o or -otar\n$usage"
	exit 1
fi
if [ "$out_manifest" = "" ]; then
	echo -e "Error: set -o or -oman\n$usage"
	exit 1
fi
if [ "$tmp_tar" = "" ]; then
	echo -e "Error: set -o or -otmptar\n$usage"
	exit 1
fi
if [ "$tmp_manifest" = "" ]; then
	echo -e "Error: set -o or -otmpman\n$usage"
	exit 1
fi
# check optional arguments
if [ "$ref_filter" = "" ]; then
	echo -e "-F not set, will not perform read filtering"
fi
# fix for Galaxy: None will be reset as empty string
if [ "$ref_filter" = "None" ]; then
	echo -e "No reference filter will be performed"
	ref_filter=""
fi
# check flex setting
if [ "$flex" != "y" ]; then
	flex=""
else
	flex="-f"
fi
# check pmethod
if [ "$pmethod" != "3p" ] && [ "$pmethod" != "5p" ] && [ "$pmethod" != "center" ]; then
	echo -e "Set -Pmethod to be '3p, '5p', or 'center'. '$pmethod' is not valid"
	exit 1
fi

# if $out is undefined, pick it based on the output manifest file
if [ "$out" = "" ]; then
	out=$(echo "$out_manifest" | sed 's/.txt//g')
fi

# get directories of the manifests

# add absolute path to out
if [[ "$out" != "/"* ]]; then
	out_dir=$(cd $(dirname $out) && pwd -P)
	base=$(basename $out)
	out="$out_dir/$base"
	out_manifest="$out.manifest.txt"
	tmp_manifest="$out.tmp.manifest.txt"
	out_tar="$out.results.tar.gz"
	tmp_tar="$out.tmp.tar.gz"
else
	out_dir=$(dirname $out)
fi

out_tar_base=$(basename $out_tar)
tmp_tar_base=$(basename $tmp_tar)

# start output manifests
echo -e "$out_tar_base\tAll\ttar.gz" | sed "s!$out_dir/!!g" > $out_manifest
echo -e "$tmp_tar_base\tAll\ttar.gz" | sed "s!$out_dir/!!g" > $tmp_manifest



#
# utility functions
#
# check if file $1 exists
# if it does, add it to the file list, and make a note in the manifest that the sample name is $3 and the type is $4
# file list and manifest based on $2
function check_file {
	if [ ! -f "$1" ]; then
	   (>&2 echo "Error: file $1 does not exist")
	   exit 1
	fi
	# append file to list
	if [ "$2" = "tmp" ]; then
		echo -e "$1\t$3\t$4" | sed "s!$out_dir/!!g" >> $tmp_manifest
	elif [ "$2" = "out" ]; then
		echo -e "$1\t$3\t$4" | sed "s!$out_dir/!!g" >> $out_manifest
	fi
}

# check on multiple files
function check_many_files {
	list=$(ls $1)
	for f in $list; do
		check_file $f $2 $3 $4
	done
}

# find executable for $1: either itself, or the $path/$1
function find {
	# check in the $path provided first
	if [ "$path" != "" ]; then
		check=$(ls $path/$1 2> /dev/null | wc -l | awk '{print $1}')
		if [ "$check" = 1 ]; then
			echo "$path/$1"
			return
		fi
	fi
	# otherwise look in user's path
	check=$(which $1 2> /dev/null | wc -l | awk '{print $1}')
	if [ "$check" = 1 ]; then
		echo "$1"
	else
		(>&2 echo "Error: script $1 can't be found")
		exit 1
	fi
}

#
# start the pipeline
#
# load dependencies
bowtie2=$(find bowtie2)
bowtie2_build=$(find bowtie2-build)
split_fastq_by_barcode_in_read=$(find split_fastq_by_barcode_in_read.pl)
hard_trim=$(find hard_trim.sh)
samtools=$(find samtools)
BAM_To_AssignCount_input_format=$(find BAM_To_AssignCount_input_format.pl)
assignCount3prime_flexible=$(find assignCount3prime_flexible.py)
GeneCount=$(find GeneCount.py)
rpm=$(find rpm.pl)
wigToBigWig=$(find wigToBigWig)

# make copy of genome
cp $ref_genome $out.genome.fa
ref_genome="$out.genome.fa"
check_file $ref_genome "tmp" "Reference" "Genome fasta"
if [ "$ref_filter" != "" ]; then
	cp $ref_filter $out.filter.fa
	ref_filter="$out.filter.fa"
	check_file $ref_filter "tmp" "Reference" "Filter fasta"
fi

# make indices
echo -e "\n#\n# Building reference indices\n#"
$bowtie2_build $ref_genome $out.genome &> /dev/null
check_many_files "$out.genome*.bt2" "tmp" "Reference" "Genome index"
if [ "$ref_filter" != "" ]; then
	$bowtie2-build $ref_filter $out.filter &> /dev/null
	check_many_files "$out.filter*.bt2" "tmp" "Reference" "Filter index"
fi
# get name of genome
$samtools faidx $ref_genome
check_file "$ref_genome.fai" "tmp" "Reference" "Fasta index"
# make 1 bigger to fit an extra base
awk -F "\t" '{OFS="\t"; $2++; print }' < $ref_genome.fai > $ref_genome.tmp && mv $ref_genome.tmp $ref_genome.fai
gname=$(head -1 $ref_genome | cut -f1 -d ' ' | sed 's/>//')

# demultiplexing
if [ -f "$barcodes" ]; then
	echo -e "\n#\n# Demultiplexing\n#"
	$split_fastq_by_barcode_in_read -i $barcodes -r1 $fastq -o $out > $out.demultiplexing_report.txt $flex
	check_file $out.demultiplexing_report.txt "out" "All" "Report"
	
	# get list of outputs
	blist=$(cat $barcodes | awk '{print $NF}')
else
	echo -e "\n#\n# NOTE: Skipping demultiplexing, all reads will belong to the same sample named 'Sample'\n#"
	blist="sample"
	ln -sf $fastq ${out}sample_R1.fastq.gz
fi

#
# processing per sample
#
# store stderr in a tmp file to echo to stdout later
echo -e "\n#\n# Per-sample processing\n#"
for b in $blist; do
	# check this file exists
	echo -e "\n#\n# Sample $b\n#"
	raw="${out}${b}_R1.fastq.gz"
	check_file $raw "tmp" $b "fastq"

	# skip if empty
	if [ -f "$barcodes" ]; then
		count=$(awk -v name="$b" '{if($1==name) print $2}' < $out.demultiplexing_report.txt)
	else
		count=$(($minreads+1))
	fi
	if (( $count >= $minreads )); then
		# start a QC file
		echo -e "Raw reads\t$count" > $out$b.stats.txt
		check_file $out$b.stats.txt "out" $b "Stats"

		# hard trim
		echo -e "# Hard trim for $b: $trim5p from 5', $trim3p from 3'"
		$hard_trim $raw $out$b.trim.fastq.gz $trim5p $trim3p
		raw="$out$b.trim.fastq.gz"
		check_file $raw "tmp" $b "fastq-trimmed"

		# rRNA filtering
		if [ "$ref_filter" != "" ]; then
			echo -e "# Filtering against filter reference $ref_filter with seed $filter_seed"
			# save alignments to get each set
			$bowtie2 -L $filter_seed -x $out.filter -U $raw 2> $out.STDERR | $samtools view -b > $out$b.filter.orig.bam
			# get unmapped reads
			$samtools view -f4 $out$b.filter.orig.bam | $samtools fastq - 2> /dev/null | gzip > $out$b.filter.fastq.gz
			# get mapped reads
			$samtools view -F4 $out$b.filter.orig.bam | $samtools fastq - 2> /dev/null | gzip > $out$b.mapped-to-filter.fastq.gz
			# count hits per filter seq
			$samtools view -F4 $out$b.filter.orig.bam | cut -f3 | sort | uniq -c | sort -gr | awk 'BEGIN{OFS="\t"; print "Filter seq","Count"}{print $2, $1}' > $out$b.mapped-to-filter.summary.txt
			rm $out$b.filter.orig.bam
			cat $out.STDERR
			raw="$out$b.filter.fastq.gz"
			check_file $raw "tmp" $b "fastq-filtered"
			check_file "$out$b.mapped-to-filter.fastq.gz" "tmp" $b "mapped-to-filter-fastq"
			check_file "$out$b.mapped-to-filter.summary.txt" "out" $b "mapped-to-filter-summary"

			# count reads
			gzip -dc $raw | wc -l | awk '{print "Reads passing alignment filter\t"$1/4}' >> $out$b.stats.txt
		fi

		# genome alignment
		echo -e "# Genome alignment to $ref_genome with seed $genome_seed"
		$bowtie2 -L $genome_seed -x $out.genome -U $raw 2> $out.STDERR | $samtools view -b - | $samtools sort - > $out$b.bam 2> /dev/null
		cat $out.STDERR
		check_file $out$b.bam "tmp" $b "BAM"
		$samtools view -F4 $out$b.bam | wc -l | awk '{print "Aligned reads\t"$1}' >> $out$b.stats.txt
	
		# assign to counts
		echo -e "# Running BAM_To_AssignCount_input_format.pl with mismatches $mismatch"
		# convert to sam to avoid samtools dependency in BAM_To_AssignCount_input_format.pl
		$samtools view -h $out$b.bam > $out$b.sam
		$BAM_To_AssignCount_input_format --input $out$b.sam --output $out$b.assignCount.txt --mismatch $mismatch
		check_file $out$b.assignCount.txt "tmp" $b "AssignCount"
		rm $out$b.sam
		wc -l $out$b.assignCount.txt | awk '{print "Reads passing mismatch filter\t"$1}' >> $out$b.stats.txt

		# assign p-site
		echo -e "# Assigning p-site based on cutoff $psite_offset, for reads from $min_len to $max_len"
		$assignCount3prime_flexible $psite_offset $out$b.assignCount.txt $out$b.posCount.txt $out$b.negCount.txt $out$b.nuclBias.txt $out$b.lengthDist.txt $min_len $max_len $pmethod
		check_file $out$b.posCount.txt "tmp" $b "Psite-posCount"
		check_file $out$b.negCount.txt "tmp" $b "Psite-negCount"
		check_file $out$b.nuclBias.txt "tmp" $b "Psite-nuclBias"
		check_file $out$b.lengthDist.txt "out" $b "Psite-lengthDist"

		# counts per gene
		echo -e "# Generating counts per gene"
		$GeneCount $out$b.posCount.txt $out$b.negCount.txt $ref_annotation $out$b.genes.txt $out$b.genesNorm.txt
		check_file $out$b.genes.txt "out" $b "Gene_counts"
		check_file $out$b.genesNorm.txt "tmp" $b "Gene_norm"

		# rpm
		echo -e "# RPM normalization"
		$rpm $out$b.genes.txt 5 $out$b.posCount.txt $out$b.negCount.txt 2 $gname $out$b.posRPM.txt $out$b.posRPM.wig $out$b.negRPM.txt $out$b.negRPM.wig
		check_file $out$b.posRPM.txt "out" $b "Gene_posRPM"
		check_file $out$b.posRPM.wig "out" $b "pos_wig"
		check_file $out$b.negRPM.txt "out" $b "Gene_negRPM"
		check_file $out$b.negRPM.wig "out" $b "neg_wig"

		# convert to bigWig
		echo -e "# Converting wig to bigWig"
		$wigToBigWig $out$b.posRPM.wig $ref_genome.fai $out$b.posRPM.bw
		$wigToBigWig $out$b.negRPM.wig $ref_genome.fai $out$b.negRPM.bw
		check_file $out$b.posRPM.bw "out" $b "pos_bw"
		check_file $out$b.negRPM.bw "out" $b "neg_bw"

	else
		echo "Skipping sample $b: insufficient data ($count reads < $minreads minimum)"
	fi
done

# pack it all up, excluding the tarball itself
# exclude the file paths when building the tarball
outlist=$(cut -f1 $out_manifest | grep -vw "$out_tar_base" | sed "s!$out_dir/!!g" | paste -s -d ' ' -)
cd $out_dir && tar -zcf $out_tar $outlist 2> $out.STDERR
cat $out.STDERR
echo $outlist | xargs -n1 rm

tmplist=$(cut -f1 $tmp_manifest | grep -vw "$tmp_tar_base" | sed "s!$out_dir/!!g" | paste -s -d ' ' -)
cd $out_dir && tar -zcf $tmp_tar $tmplist 2> $out.STDERR
cat $out.STDERR
echo $tmplist | xargs -n1 rm

# remove the STDERR report
rm $out.STDERR

exit 0
