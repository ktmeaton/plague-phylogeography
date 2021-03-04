#!/bin/bash

RESULTS_DIR=$1
MODE=$2

KEEP_DIR=(
	data
	sqlite_db
	eager
	snippy_pairwise
	qualimap
	);

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
	if [[ $keep == "false" && $MODE == "list" ]];
	then
			echo -e "\tRemoving: ${RESULTS_DIR}/$dirname";
			echo -e "\t         rm -rf ${RESULTS_DIR}/$dirname";
	elif [[ $keep == "false" && $MODE == "rm" ]];
	then
			echo -e "\tRemoving: ${RESULTS_DIR}/$dirname";
			echo -e "\t         rm -rf ${RESULTS_DIR}/$dirname";
			rm -rf ${RESULTS_DIR}/$dirname
	fi;
done;
