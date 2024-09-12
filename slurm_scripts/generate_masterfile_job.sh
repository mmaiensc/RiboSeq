#!/bin/bash
#SBATCH -J masterfile
#SBATCH -t 100:00:00
#SBATCH -N 1
#SBATCH -n 1
######################################
# EDIT BELOW HERE TO CHANGE SETTINGS #
######################################
# for help, see:
# https://github.com/mmaiensc/RiboSeq/wiki/Generating-master-files
DATADIR="/projects/psci_shura_chi/riboseq/data"
# required
CONTROL_MANIFEST="riboseq_practice.manifest.txt"
TREATMENT_MANIFEST="riboseq_practice.manifest.txt"
CONTROL_SAMPLE="Sample1"
TREATMENT_SAMPLE="Sample2"
REFERENCE="practice_reference.fa"
ANNOTATION="practice_annotation.txt"
OUTPUT_NAME="riboseq_practice_sample1-sample2.masterfile.txt"

# optional
## set filter to emtpy (FILTER="") to skip
GENE_THRESHOLD="100"
CODON_THRESHOLD_SAMPLE="5"
CODON_THRESHOLD_TREATMENT="5"
RPM_SEGMENT="segment"

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
SCRIPTDIR="/projects/psci_shura_chi/riboseq/RiboSeq"

module purge
module load SAMtools &> /dev/null

# set up filter
$SCRIPTDIR/generate_master.sh \
	-m1 $DATADIR/$CONTROL_MANIFEST \
	-i1 $CONTROL_SAMPLE \
	-m2 $DATADIR/$TREATMENT_MANIFEST \
	-i2 $TREATMENT_SAMPLE \
	-G $DATADIR/$REFERENCE \
	-A $DATADIR/$ANNOTATION \
	-o $DATADIR/$OUTPUT_NAME \
	-c $GENE_THRESHOLD \
	-c1 $CODON_THRESHOLD_SAMPLE \
	-c2 $CODON_THRESHOLD_TREATMENT \
	-norm $RPM_SEGMENT \
	-path $SCRIPTDIR &> $DATADIR/$OUTPUT_NAME.log

