#!/bin/bash

INFILE=$1
QCHAR="I"

awk -v QCHAR=$QCHAR 'BEGIN{
  RS=">";
	ORS="\n";
	FS="\n";
	i=1;
}
{
  if (NR > 1)
	{
		print "@seq"i;
		print $2;
		print "+";
		qual = ""
		for(c=0;c<length($2);c++){ qual = qual QCHAR }
		print qual;
		i+=1
  }
}' $INFILE
