#!/bin/bash

RESULTS_DIR=$1
CONFIRM=$2

KEEP_DIR=(data sqlite_db);

if [[ ! $RESULTS_DIR ]];
then
	exit 1;
fi

echo -e "\nCleaning ${RESULTS_DIR} of the following directories: "

for dir in `ls -d $RESULTS_DIR/*`;
do
	dirname=`basename $dir`;
	# Iterate through target directories to keep
	keep="false";
	for target in ${KEEP_DIR[@]};
	do
		if [[ $dirname == $target ]];
		then
			keep="true";
	  fi;
	done;

	# Check if this directory should be kept
	if [[ $keep == "false" ]];
	then
		echo $dirname
	fi
done;
