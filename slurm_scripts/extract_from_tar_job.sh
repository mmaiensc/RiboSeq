#!/bin/bash
#SBATCH -J extract_from_tar
#SBATCH -t 100:00:00
#SBATCH -N 1
#SBATCH -n 1
######################################
# EDIT BELOW HERE TO CHANGE SETTINGS #
######################################
# for help, see:
# https://github.com/mmaiensc/RiboSeq/wiki/Extracting-files-from-tar-archives
DATADIR="/projects/psci_shura/riboseq/data"
# required
TAR="riboseq_practice.results.tar.gz"
WANTED_FILE="riboseq_practiceSample1.stats.txt"
OUTPUT_NAME="riboseq_practice_sample_stats.txt"

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

# run command
$SCRIPTDIR/extract_from_tar.sh $DATADIR/$TAR $WANTED_FILE $DATADIR/$OUTPUT_NAME &> $DATADIR/$OUTPUT_NAME.log

