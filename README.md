# paper-phylogeography
Phylogeography of Yersinia pestis

## Dependencies
**Workflow:** NextFlow  
**Database:** NCBImeta, sqlite3 (CLI)  
**Alignment:** snippy  
**Masking, etc.:** dustmasker, mummer  

### Conda Environment
Create a conda environment with the required dependencies  
```
conda env create -f phylo-env.yaml --name phylo-env
conda activate phylo-env
```

### Full Pipeline
```
nextflow run pipeline.nf \
  --ncbimeta_create ncbimeta.yaml \
  --ncbimeta_update ncbimeta.yaml \
  -with-trace
  -with-timeline
  -with-dag pipeline.pdf
  -with-report
```

## Fresh start

### Build NCBImeta database, test lite run-through of pipeline
```
nextflow run pipeline.nf \
  --ncbimeta_create ncbimeta.yaml \
  --skip_sqlite_import
```

### Annotate the Database
Query the Database for problematic records (wrong organism)
```
DB=results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite
sqlite3 $DB
.output extract.txt
SELECT BioSampleAccession,
       BioSampleBioProjectAccession,
       BioSampleStrain,
       BioSampleOrganism,
       BioSampleSRAAccession,
       BioSampleAccessionSecondary,
       BioSampleCollectionDate,
       BioSampleGeographicLocation,
       BioSampleHost,
       BioSampleComment
FROM BioSample
WHERE (BioSampleOrganism NOT LIKE '%Yersinia pestis%');
```
Add delimited headers to top of file
```
DELIM="|"
sed  -i "1i BioSampleAccession${DELIM}BioSampleBioProjectAccession${DELIM}BioSampleStrain${DELIM}BioSampleOrganism${DELIM}BioSampleSRAAccession${DELIM}BioSampleAccessionSecondary${DELIM}BioSampleCollectionDate${DELIM}BioSampleGeographicLocation${DELIM}BioSampleHost${DELIM}BioSampleComment" extract.txt
```
Convert from pipe-separated to tab-separated file
```
sed -i "s/|/\t/g" extract.txt
```
Inspect the extract.txt file in a spreadsheet view (ex. Excel, Google Sheets)  
Add "REMOVE: Not Yersinia pestis" to the BioSampleComment column to any rows that are confirmed appropriate.  


### Update Database With Annotations
```
nextflow run pipeline.nf \
  --ncbimeta_update ncbimeta.yaml \
  --ncbimeta_annot extract.txt \
  --skip_sqlite_import \
  -resume
```

### Run from established database
```
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --ncbimeta_annot annot.txt \
  --max_datasets 1 \
  -resume
```

### Tracing
-with-trace
-with-timeline
-with-dag pipeline.pdf
-with-report
