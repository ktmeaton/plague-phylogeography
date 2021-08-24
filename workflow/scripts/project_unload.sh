#!/bin/bash

RESULTS_DIR=$1
BACKUP_DIR=$2
MODE=$3

EXCLUDE_DIR=(
	# Input dir
	#config
	data
	sqlite_db
	# Output per-sample
	detect_snp_density
	dnds
	eager
	locus_coverage
	qualimap
	snippy_pairwise
	);

if [[ ! $RESULTS_DIR || ! $BACKUP_DIR ]];
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

    # Print directory specific info
	if [[ $keep == "true" ]]; then
        echo -e "\tBacking up: ${RESULTS_DIR}/$dirname/";
        echo -e "\t    origin: ${RESULTS_DIR}/$dirname/";
        echo -e "\t      dest: ${BACKUP_DIR}/$dirname/";
    else
        continue
    fi

    # Print MODE specific output
	if [[ $MODE == "cp" ]]; then
		echo -e "\t         cp -r ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname";
		cp -r ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname
	elif [[ $MODE == "mv" ]]; then
		echo -e "\t         mv ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname";
		mv ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname
	elif [[ $MODE == "rsync" ]]; then
		echo -e "\t       rsync -u -a ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/$dirname/";
		rsync -u -a ${RESULTS_DIR}/$dirname ${BACKUP_DIR}/;
	fi
    echo

done;
