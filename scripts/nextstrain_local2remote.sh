#/bin/bash

# nextstrain_local2remote.sh User Repository Narrative Dataset
# nextstrain_local2remote.sh ktmeaton plague-phylogeography plagueSCDS2020 plague150

USER=$1
REPO=$2
NARRATIVE=$3
DATA=$4

# Defaults
LOCAL_URL="http://localhost:4000"
REMOTE_URL="https://nextstrain.org/community"
NARRATIVE_DIR="narratives"


# First change data path to remote, then change narrative path to remote
sed "s|${LOCAL_URL}/${DATA}Local|${REMOTE_URL}/${USER}/${REPO}/${DATA}Remote|g" ${NARRATIVE_DIR}/${NARRATIVE}Local.md | \
	sed "s|${LOCAL_URL}/${NARRATIVE_DIR}/${NARRATIVE}Local|${REMOTE_URL}/${NARRATIVE_DIR}/${USER}/${REPO}/${NARRATIVE}Remote|g" \
	> ${NARRATIVE_DIR}/${REPO}_${NARRATIVE}Remote.md
