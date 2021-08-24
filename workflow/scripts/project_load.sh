#!/bin/bash

RESULTS_DIR=$1
BACKUP_DIR=$2
MODE=$3

# Credits: @Dave Dopson
# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# All directories
EXCLUDE_DIR=(
  #beast
	);

if [[ ! $RESULTS_DIR || ! $BACKUP_DIR ]];
then
	echo "Please specify a results directory and backup directory."
	exit 1;
fi

echo -e "Cleaning directory ${RESULTS_DIR} before loading..."
echo
#${SCRIPT_DIR}/project_clean.sh ${RESULTS_DIR}

echo -e "\nLoading project ${BACKUP_DIR} to ${RESULTS_DIR}:"
echo

mkdir -p ${RESULTS_DIR}

for dir in `ls -d $BACKUP_DIR/*`;
do
    if [[ ! -d $dir ]]; then
		continue
	fi
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
        echo -e "\tLoading: ${BACKUP_DIR}/$dirname/";
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
