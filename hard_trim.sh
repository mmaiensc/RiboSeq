# hard trim N bases from start and N bases from end of read
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] || [ "$4" = "" ]; then
	echo "Usage: input.fastq.gz output.fastq.gz N_start N_end"
	echo "  Trims N_start from start of read, N_end from end of read"
	echo "  Only prints reads if length is >5 after trimming"
	exit
fi

# run hard trimming
gzip -dc $1 | paste - - - - | awk -v s="$3" -v e="$4" -F "\t" '{if(length($2) > s+e+5){$2=substr($2, s+1, length($2)-e-s); $4=substr($4, s+1, length($4)-e-s); print $1; print $2; print $3; print $4}}' | gzip > $2



