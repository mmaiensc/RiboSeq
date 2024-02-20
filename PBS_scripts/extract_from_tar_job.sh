#PBS -N extract_from_tar
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
######################################
# EDIT BELOW HERE TO CHANGE SETTINGS #
######################################
# for help, see:
# https://github.com/mmaiensc/RiboSeq/wiki/Extracting-files-from-tar-archives
DATADIR="/mnt/store2/home/mmaiensc/ondemand/ribo_profiling/data"
# required
TAR="practice.results.tar.gz"
WANTED_FILE="practiceSample1.stats.txt"
OUTPUT_NAME="practice_sample_stats.txt"

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

# run command
$SCRIPTDIR/extract_from_tar.sh $DATADIR/$TAR $WANTED_FILE $DATADIR/$OUTPUT_NAME &> $DATADIR/$OUTPUT_NAME.log

