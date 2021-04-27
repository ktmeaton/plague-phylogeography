#!/bin/bash

# File variables
INFILE=$1

# Create bed format using awk processing
awk 'BEGIN{RS = ">"; FS = "\n"}
{
	if(NR > 1)
	{
		fasta_seq = ""
		for (i=2; i<=NF; i++)
		{
			fasta_seq = fasta_seq$i
		}
		print ">"$1"\n"fasta_seq;

	}

}' $INFILE
