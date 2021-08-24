#!/bin/bash

FILE=$1

awk '{ \
if (NR == 3) \
{ \
    split($0,sampleArr,","); \
} \
if (NR == 6) \
{ \
    split($0,splitArr,","); \
} \
} \
END{ \
for (i=2; i<= length(sampleArr); i++) \
{ \
print sampleArr[i] "\t" splitArr[i]; \
}}' $FILE
