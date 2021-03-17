#!/bin/bash

RESULTS_DIR=$1
BACKUP_DIR=$2
MODE=$3

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

if [[ ! $RESULTS_DIR || ! $BACKUP_DIR ]];
then
	exit 1;
fi

echo -e "\nRestoring ${BACKUP_DIR} to ${RESULTS_DIR}:"
echo
mkdir -p ${RESULTS_DIR}

for dir in `ls -d $BACKUP_DIR/*`;
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
        echo -e "\tRestoring: ${BACKUP_DIR}/$dirname/";
        echo -e "\t   origin: ${BACKUP_DIR}/$dirname/";
        echo -e "\t     dest: ${RESULTS_DIR}/$dirname/";
    else
        continue
    fi

    # Print MODE specific output
	if [[ $MODE == "cp" ]]; then
		echo -e "\t         cp -r ${BACKUP_DIR}/$dirname ${RESULTS_DIR}/"
        cp -r ${BACKUP_DIR}/$dirname ${RESULTS_DIR}/;
	elif [[ $MODE == "mv" ]]; then
		echo -e "\t        mv ${BACKUP_DIR}/$dirname ${RESULTS_DIR}/"
        mv ${BACKUP_DIR}/$dirname ${RESULTS_DIR}/;
	elif [[ $MODE == "rsync" ]]; then
		echo -e "\t           rsync -u -a ${BACKUP_DIR}/$dirname ${RESULTS_DIR}/ ";
        rsync -u -a ${BACKUP_DIR}/$dirname ${RESULTS_DIR}/;
	fi
    echo
done;
