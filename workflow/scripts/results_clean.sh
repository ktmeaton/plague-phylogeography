#!/bin/bash

RESULTS_DIR=$1
MODE=$2

EXCLUDE_DIR=(
    config
	data
	detect_snp_density
	eager
	locus_coverage
    logs
	qualimap
	snippy_pairwise
	sqlite_db
	);

if [[ ! $RESULTS_DIR ]];
then
	exit 1;
fi

echo -e "\nCleaning up ${RESULTS_DIR} directory:"
mkdir -p $BACKUP_DIR;

for dir in `ls -d $RESULTS_DIR/*`;
do
	dirname=`basename $dir`;
	# Iterate through target directories to exclude
	keep="true";
	for target in ${EXCLUDE_DIR[@]};
	do
		if [[ $dirname == $target ]]; then
			keep="false";
	    fi;
	done;

    # Print directory specific info
	if [[ $keep == "true" ]]; then
        echo -e "\tCleaning up: ${RESULTS_DIR}/$dirname/";
        echo -e "\t     origin: ${RESULTS_DIR}/$dirname/";
        echo -e "\t       dest: ${RESULTS_DIR}/$dirname/";
    else
        continue
    fi

	# Delete Directories
	if [[ $MODE == "rm" ]]; then
		echo -e "\t             rm -rf ${RESULTS_DIR}/$dirname/";
		rm -rf ${RESULTS_DIR}/$dirname/;
	fi

	echo

done;
