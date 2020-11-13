#!/bin/bash

# Output Format
# 1. BioSampleAccession
# 2. Strain
# 3. Country
# 4. Province
# 6. BeginDate
# 7. EndDate
# 8. LatLon

QUERY="SELECT
  BioSampleAccession,
	BioSampleStrain,
  BioSampleGeographicLocation,
  BioSampleCollectionDate
FROM
  BioSample
WHERE
  BioSampleComment LIKE '%KEEP%SRA%Ancient'";

DELIM="\t"

echo -e "Accession"$DELIM"Strain"$DELIM"Country"$DELIM"Province"$DELIM"Begin"$DELIM"End"$DELIM"LatLon"

sqlite3 \
  results/sqlite_db/yersinia_pestis_db.sqlite \
  "$QUERY" | \
  while read line;
  do
    acc=`echo $line | cut -d "|" -f 1`;
		strain=`echo $line | cut -d "|" -f 2`;
    country=`echo $line | cut -d "|" -f 3 | cut -d ":" -f 1`;
    province=`echo $line | cut -d "|" -f 3 | cut -d ":" -f 2`;
    begin=`echo $line | cut -d "|" -f 4 | cut -d ":" -f 1`;
    end=`echo $line | cut -d "|" -f 4 | cut -d ":" -f 2`;
		latlon=`python workflow/scripts/geocode.py "$country:$province"`;
		echo -e	$acc$DELIM$strain$DELIM$country$DELIM$province$DELIM$begin$DELIM$end$DELIM$latlon;
  done
