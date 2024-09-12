#!/bin/bash
#SBATCH -J ribo_profiling
#SBATCH -t 100:00:00
#SBATCH -N 1
#SBATCH -n 1
######################################
# EDIT BELOW HERE TO CHANGE SETTINGS #
######################################
# for help, see:
# https://github.com/mmaiensc/RiboSeq/wiki/RiboSeq-profiling
DATADIR="/projects/psci_shura/riboseq/data"
# required
DATA="practice.fastq.gz"
BARCODES="practice_barcodes.txt"
REFERENCE="practice_reference.fa"
ANNOTATION="practice_annotation.txt"
OUTPUT_NAME="riboseq_practice"

# optional
## set filter to emtpy (FILTER="") to skip
FILTER="practice_filter.fa"
min="1000"
flex="n"
hard5="2"
hard3="5"
Fseed="10"
Gseed="18"
mis="2"
Poff="15"
Pmin="24"
Pmax="46"

################
# STOP EDITING #
################


#
# Move to the directory where the job was submitted
if [ "$SLURM_SUBMIT_DIR" != "" ]; then
	cd $SLURM_SUBMIT_DIR
fi
#
echo "Job was submitted from: " $SLURM_SUBMIT_DIR

# paths to scripts, data, outputs
SCRIPTDIR="/projects/psci_shura/riboseq/RiboSeq"

# modules
module purge
module load Bowtie2
module load BEDTools &> /dev/null
module load SAMtools &> /dev/null
module load kent_tools

# set up filter
if [ "$FILTER" != "" ]; then
	FILTER="-F $DATADIR/$FILTER"
fi

$SCRIPTDIR/ribo-seq_sample_processing.sh \
	-f $DATADIR/$DATA \
	-b $DATADIR/$BARCODES \
	-G $DATADIR/$REFERENCE \
	-A $DATADIR/$ANNOTATION \
	-o $DATADIR/$OUTPUT_NAME \
	-min $min \
	-flex $flex \
	-hard5 $hard5 \
	-hard3 $hard3 \
	-Fseed $Fseed \
	-Gseed $Gseed \
	-mis $mis \
	-Poff $Poff \
	-Pmin $Pmin \
	-Pmax $Pmax \
	$FILTER -path $SCRIPTDIR &> $DATADIR/$OUTPUT_NAME.log

