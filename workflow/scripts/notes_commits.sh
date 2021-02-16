#!/bin/bash

# First ="a58299c51609ad"
OLD_VER=$1
NEW_VER=$2
BASE="https://github.com/ktmeaton/ActionRPG/commit"

git log --pretty=oneline --abbrev-commit ${OLD_VER}..${NEW_VER} | \
  while read line;
  do
    hash=`echo $line | cut -d " " -f 1`
    msg=`echo $line | sed "s/$hash //g"`
    echo -e "* [\`\`\`$hash\`\`\`]($BASE/$hash) $msg";
  done;
