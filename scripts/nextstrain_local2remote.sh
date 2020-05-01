#/bin/bash

#MARKDOWN_LOCAL="narratives/plagueSCDS2020Local.md";
#MARKDOWN_REMOTE="narratives/plagueSCDS2020Remote.md";
#DATA_LOCAL="http://localhost:4000/plague150Local";
#DATA_REMOTE="https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote";
#NARRATIVE_LOCAL="http://localhost:4000/narratives/plagueSCDS2020Local";
#NARRATIVE_REMOTE="https://nextstrain.org/community/narratives/ktmeaton/plague-phylogeography/plagueSCDS2020Remote";

# nextstrain_local2remote.sh ktmeaton/paper-phylogeography
#s|http://localhost:4000/plague150Local|https://nextstrain.org/community/ktmeaton/plague-phylogeography/plague150Remote|g

USER="ktmeaton"
REPO="plague-phylogeography"
NARRATIVE="plagueSCDS2020"
DATA="plague150"

LOCAL_URL="http://localhost:4000"
REMOTE_URL="https://nextstrain.org/community"
NARRATIVE_DIR="narratives"

echo ${LOCAL_URL}${NARRATIVE_DIR}${NARRATIVE}Local

# First change data path to remote, then change narrative path to remote
sed "s|${LOCAL_URL}/${DATA}Local|${REMOTE_URL}/${USER}/${REPO}/${DATA}Remote|g" ${NARRATIVE_DIR}/${NARRATIVE}Local.md | \
	sed "s|${LOCAL_URL}/${NARRATIVE_DIR}/${NARRATIVE}Local|${REMOTE_URL}/${NARRATIVE_DIR}/${USER}/${REPO}/${NARRATIVE}Remote|g" \
	> ${NARRATIVE_DIR}/${REPO}_${NARRATIVE}Remote.md


# sed "s|${DATA_LOCAL}|${DATA_REMOTE}|g" ${MARKDOWN_LOCAL} | sed "s|${NARRATIVE_LOCAL}|${NARRATIVE_REMOTE}|g" 


