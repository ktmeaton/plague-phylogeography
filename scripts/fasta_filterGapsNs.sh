#!/bin/bash

# Assumes unwrapped fasta, no spaces between records
# Assumes missing/ambig data represented as "-","N", or "X"

# File variables
INFILE=$1
MISSINGFRAC=$2
BACKBONEFILE=$3

# Remove duplicates using awk processing
awk -v missFrac="$MISSINGFRAC" -v backBoneFile="$BACKBONEFILE" '{
     if(NR % 2 == 1)
     {
         header=$0
     }
     if(NR == 2)
     {
         for(i=1; i<= length($0); i++)
         {
             excludeArr[i] = 0;
         }
     }
     if(NR % 2 == 0)
     {
         seqArr[header]=$0;
         split($0, seq, "");
         for (i=1; i <= length($0); i++)
         {
                 if(seq[i] == "-" || toupper(seq[i]) == "N" || toupper(seq[i]) == "X")
                 {
                     excludeArr[i]+= 1;
                 }
         }
     }
     } END{
     firstSeq= 1;
     for (header in seqArr)
     {
         print header;
         split(seqArr[header], seq, "");
         finalSeq = "";
         for (i=1; i <= length(seqArr[header]); i++)
         {
              if (excludeArr[i]  / (NR / 2) <= missFrac)
              {
                  finalSeq = finalSeq seq[i];
                  if(firstSeq == 1)
                  {
                      if (!beginBackBone)
                      {
                           beginBackBone=i;
                           endBackBone=i;
                      }
                      else if (i - 1 == endBackBone)
                      {
                          endBackBone=i;
                      }
                      else if (i - 1 > endBackBone)
                      {
                          print beginBackBone"-"endBackBone >> backBoneFile;
                          beginBackBone=i;
                          endBackBone=i;
                      }
                  }
              }
         }
         if(firstSeq == 1){print beginBackBone"-"endBackBone >> backBoneFile;}
         print finalSeq;
         firstSeq = 0;
     }
     }' $INFILE
