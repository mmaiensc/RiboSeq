#PBS -N ribo_profiling
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
######################################
# EDIT BELOW HERE TO CHANGE SETTINGS #
######################################
# for help, see:
# https://github.com/mmaiensc/RiboSeq/wiki/RiboSeq-profiling
DATADIR="/mnt/store2/home/mmaiensc/ondemand/ribo_profiling/data"
# required
DATA="practice.fastq.gz"
BARCODES="practice_barcodes.txt"
REFERENCE="practice_reference.fa"
ANNOTATION="practice_annotation.txt"
OUTPUT_NAME="practice"

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
if [ "$PBS_O_WORKDIR" != "" ]; then
	cd $PBS_O_WORKDIR
fi
#
echo "Job was submitted from: " $PBS_O_WORKDIR

# paths to scripts, data, outputs
SCRIPTDIR="/mnt/store2/home/mmaiensc/ondemand/ribo_profiling/RiboSeq"

# modules -- for SABER; commented-out pertain to extreme
module purge
#module load Bowtie2
module load apps/bowtie2-2.3.5
#module load BEDTools &> /dev/null
module load SAMtools &> /dev/null
#module load cri-centos7/ucsc_tools
module load apps/ucsc_tools

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

