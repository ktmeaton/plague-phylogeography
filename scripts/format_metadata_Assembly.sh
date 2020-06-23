#!/bin/bash

# 1st Param: Pipeline outdir
#    Ex. Assembly_Modern_Outgroup
# 2nd Param: SQLite DB Path
#    Ex. ~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite
# 3rd Param: Project script directory
#    Ex. ~/.nextflow/assets/ktmeaton/plague-phylogeography/scripts/
#
# Example Usage:
# format_metadata_Assembly.sh \
#   Assembly_Modern_Outgroup \
#   ~/.nextflow/assets/ktmeaton/plague-phylogeography/results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
#   ~/.nextflow/assets/ktmeaton/plague-phylogeography/scripts/

# Extract metadata from sqlite database
mkdir -p $project/nextstrain/
$scriptsDir/sqlite_NextStrain_tsv.py \
  --database $sqliteDB \
  --query "SELECT BioSampleAccession,AssemblyFTPGenbank,BioSampleStrain,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleBiovar,BioSampleHost FROM Master WHERE (BioSampleComment LIKE '%KEEP%Assembly%')" \
  --no-data-char ? \
  --output $project/nextstrain/metadata_nextstrain.tsv

# Break up semi-colon separated FTP (two versions of assembly)
awk -F "\t" -v colFTP=2 'BEGIN{OFS=FS}{
if ($colFTP ~ /;/){
  split($colFTP,splitFTP,";");
  for (i=1; i<=length(splitFTP); i++){
    $colFTP=splitFTP[i];
    print $0
  }
}
else{
  print $0;
}}' $project/nextstrain/metadata_nextstrain.tsv > $project/nextstrain/metadata_nextstrain_split.tsv


# Add the reference genome metadata as a final line
sqlite3 $sqliteDB  \
  "SELECT BioSampleAccession,AssemblyFTPGenbank,BioSampleStrain,BioSampleCollectionDate,BioSampleGeographicLocation,BioSampleBiovar,BioSampleHost FROM Master WHERE BioSampleComment LIKE '%Reference%'" | \
  sed 's/|/\t/g' >> $project/nextstrain/metadata_nextstrain_split.tsv

# Write header to a new edited metadata file, add col "strain"
head -n 1 $project/nextstrain/metadata_nextstrain_split.tsv | \
  awk -F "\t" '{print "strain\t"$0}' \
  > $project/nextstrain/metadata_nextstrain_split_strain.tsv

# Figure out the assembly file names by parsing the FTP url column, save to col "strain"
tail -n +2 $project/nextstrain/metadata_nextstrain_split.tsv  | \
  awk -F "\t" '{split($2,ftpSplit,"/"); name=ftpSplit[10]"_genomic"; print name"\t"$0}' \
  >> $project/nextstrain/metadata_nextstrain_split_strain.tsv

# Change reference genome file name to "Reference"
sed -i 's/GCA_000009065.1_ASM906v1_genomic/Reference/g' $project/nextstrain/metadata_nextstrain_split_strain.tsv
# Standardize biovar nomenclature
sed -i 's/Mediaevalis/Medievalis/g' $project/nextstrain/metadata_nextstrain_split_strain.tsv
