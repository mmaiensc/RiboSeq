#PBS -N masterfile
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
######################################
# EDIT BELOW HERE TO CHANGE SETTINGS #
######################################
# for help, see:
# https://github.com/mmaiensc/RiboSeq/wiki/Generating-master-files
DATADIR="/mnt/store2/home/mmaiensc/ondemand/ribo_profiling/data"
# required
CONTROL_MANIFEST="practice.manifest.txt"
TREATMENT_MANIFEST="practice.manifest.txt"
CONTROL_SAMPLE="Sample1"
TREATMENT_SAMPLE="Sample2"
REFERENCE="practice_reference.fa"
ANNOTATION="practice_annotation.txt"
OUTPUT_NAME="practice_sample1-sample2.masterfile.txt"

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
if [ "$PBS_O_WORKDIR" != "" ]; then
	cd $PBS_O_WORKDIR
fi
#
echo "Job was submitted from: " $PBS_O_WORKDIR

# paths to scripts, data, outputs
SCRIPTDIR="/mnt/store2/home/mmaiensc/ondemand/ribo_profiling/RiboSeq"
DATADIR="/mnt/store2/home/mmaiensc/ondemand/ribo_profiling/outputs"

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

