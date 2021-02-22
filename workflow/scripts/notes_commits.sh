#!/bin/bash

# First ="de952505c2a"
OLD_VER=$1
NEW_VER=$2
BASE="https://github.com/ktmeaton/plague-phylogeography/commit"

# Default number of commits to print is 20
# Use -1 to print infinite
MAX_COMMITS=$3
MAX_COMMITS="${MAX_COMMITS:=20}"

i=0;

git log --pretty=oneline --abbrev-commit ${OLD_VER}..${NEW_VER} | \
  while read line;
  do
    if [[ $i -ge $MAX_COMMITS && $MAX_COMMITS != -1 ]]; then
      break
    fi
    hash=`echo $line | cut -d " " -f 1`
    msg=`echo $line | sed "s/$hash //g"`
    echo -e "* [\`\`\`$hash\`\`\`]($BASE/$hash) $msg";

    # Increment commit counter
    i=`expr $i + 1`
  done;
