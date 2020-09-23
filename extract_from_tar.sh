#!/bin/bash
# extract a specified file from a tar archive

if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ]; then
	echo "Usage: input.tar(.gz) file_name output"
	exit
fi

# check if we'll made a folder when extracting $2
rmpath=0
if [[ "$2" == *"/"* ]]; then
	path=$(echo $2 | awk -F "/" '{$NF=""; OFS="/"; print}')
	if [ ! -d "$path" ]; then
		rmpath=1
	fi
fi

# extract file
tar -xf $1 $2

# rename file
mv $2 $3

# delete the folder if we made it
if [ "$rmpath" = 1 ] && [ -d "$path" ]; then
	rm -r $path
fi


