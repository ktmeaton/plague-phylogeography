#!/bin/bash

FASTA=$1
OUTDIR=$2
THRESHOLD=$3
LENGTH=$4


INPREFIX=${FASTA%.*}
OUTPREFIX="${OUTDIR}/`basename $INPREFIX`"

#------------------------------------------------------------------------------#
# Inexact repeats
#------------------------------------------------------------------------------#

# Align reference to itself to find inexact repeats
nucmer --maxmatch --nosimplify --prefix=${OUTPREFIX}.inexact ${FASTA} ${FASTA}
# Convert the delta file to a simplified, tab-delimited coordinate file
show-coords -r -c -l -T ${OUTPREFIX}.inexact.delta | tail -n+5 > ${OUTPREFIX}.inexact.coords
# Remove all "repeats" that are simply each reference aligned to itself
# also retain only repeats with more than 90% sequence similarity.
awk -F "\t" -v threshold=$THRESHOLD '{if ($1 == $3 && $2 == $4 && $12 == $13)
      {next;}
  else if ($7 >= threshold)
      {print $0}}' ${OUTPREFIX}.inexact.coords > ${OUTPREFIX}.inexact.repeats
# Convert to bed file format, changing to 0-base position coordinates
awk -F "\t" '{print $12 "\t" $1-1 "\t" $2-1;
  if ($3 > $4){tmp=$4; $4=$3; $3=tmp;}
  print $13 "\t" $3-1 "\t" $4-1;}' ${OUTPREFIX}.inexact.repeats | \
sort -k1,1 -k2,2n | \
bedtools merge > ${OUTPREFIX}.inexact.repeats.bed

#------------------------------------------------------------------------------#
# Exact repeats
#------------------------------------------------------------------------------#
#echo "repeat-match"
#repeat-match -n ${LENGTH} ${FASTA} > ${OUTPREFIX}.exact.repeats
# Convert to bed file format, changing to 0-base position coordinates

#------------------------------------------------------------------------------#
# Tandem repeats
#------------------------------------------------------------------------------#
#echo "exact-tandems"
#repeat-match -t -n ${LENGTH} ${FASTA} | \
#  tail +3 | \
#  sort -k1n -k2n | \
#  awk 'BEGIN{printf "%8s %8s %8s %10s\n", "Start", "Extent", "UnitLen", "Copies";}
#        {
#         if  ($1 + $3 < $2)
#             next;
#         if  ($1 == prev)
#             next;
#         start = $1;
#         extent = $2 + $3 - $1;
#         unitlen = $2 - $1;
#         printf "%8d %8d %8d %10.1f\n", start, extent, unitlen, extent / unitlen;
#         prev = $1;
#       }' | \
#  tail -n+2 | \
#  awk -F "\t"> ${OUTPREFIX}.tandem.repeats
# Convert to bed file format, changing to 0-base position coordinates
