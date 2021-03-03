#!/bin/bash

RESULTS_DIR=$1
BACKUP_DIR=$2
MODE=$3

EXCLUDE_DIR=(
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

echo -e "\nBacking up ${RESULTS_DIR} to ${BACKUP_DIR}:"
mkdir -p $BACKUP_DIR;

for dir in `ls -d $RESULTS_DIR/*`;
do
	dirname=`basename $dir`;
	# Iterate through target directories to exclude
	keep="true";
	for target in ${EXCLUDE_DIR[@]};
	do
		if [[ $dirname == $target ]];
		then
			keep="false";
	  fi;
	done;
	# Check if this directory should be kept
	if [[ $keep == "true" && $MODE == "cp" ]];
	then
		echo -e "\tCopying: ${RESULTS_DIR}/$dirname";
		echo -e "\t         cp -r ${RESULTS_DIR}/$dirname to ${BACKUP_DIR}/$dirname";
		cp -r ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname
	elif [[ $keep == "true" && $MODE == "mv" ]];
	then
		echo -e "\tMoving: ${RESULTS_DIR}/$dirname";
		echo -e "\t         mv ${RESULTS_DIR}/$dirname to ${BACKUP_DIR}/$dirname";
		mv ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname
	elif [[ $keep == "true" && $MODE == "list" ]];
	then
		echo -e "\tMoving/Copying: ${RESULTS_DIR}/$dirname";
		echo -e "\t         mv/cp ${RESULTS_DIR}/$dirname to ${BACKUP_DIR}/$dirname";
	fi

done;
