#!/bin/bash

# Convert dustmasker interval output to bed format
# Convert 1-base coordinates to 0-base coordinations

#  interval2bed.sh Ypestis.dustmasker.intervals Ypestis.dustmasker.bed

# INFILE is the dustmasker.intervals file (dustmasker -outfmt interval)
INFILE=$1
# OUTFILE is the output bed file (dustmasker.bed)
OUTFILE=$2

awk -F "\t" 'BEGIN{RS=">";FS="\n"}
    {if (NR < 2)
        {next;}
    split($1,reference," "); 
    for (i=2; i<=NF-1;i++)
    {
        split($i,splitCoord," - ");
        startCoord=splitCoord[1]-1;
        endCoord=splitCoord[2]-1;
        print reference[1] "\t" startCoord "\t" endCoord;
    }
 }' $INFILE | sort -k1,1 -k2,2n > $OUTFILE
