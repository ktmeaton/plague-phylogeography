# paper-phylogeography
Phylogeography of Yersinia pestis

## Dependencies
**Workflow:** NextFlow  
**Database:** NCBImeta, sqlite3 (CLI)  
**Alignment:** snippy  
**Masking, etc.:** dustmasker, mummer, vcftools   
**Model Selection:** modeltest-ng  

### Conda Environment
Create a conda environment with the required dependencies  
```
conda env create -f phylo-env.yaml --name phylo-env
conda activate phylo-env
```

## Full Pipeline (Reproduce)
```
nextflow run pipeline.nf \
  --ncbimeta_create ncbimeta.yaml \
  --ncbimeta_update ncbimeta.yaml \
  --ncbimeta_annot annot_biosample.txt \
  --max_datasets 3 \
  -with-trace \
  -with-timeline \
  -with-dag pipeline.pdf \
  -with-report
```

## Step By Step (From Scratch)

### Build NCBImeta database
```
nextflow run pipeline.nf \
  --ncbimeta_create ncbimeta.yaml
```

### Annotate the Database
Query the Database for problematic records (wrong organism)
```
DB=results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite
sqlite3 $DB
.output annot_biosample.txt
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
Add delimited headers to top of file (that match NCBImeta table BioSample)
```
DELIM="|";
sed  -i "1i BioSampleAccession${DELIM}BioSampleBioProjectAccession${DELIM}BioSampleStrain${DELIM}BioSampleOrganism${DELIM}BioSampleSRAAccession${DELIM}BioSampleAccessionSecondary${DELIM}BioSampleCollectionDate${DELIM}BioSampleGeographicLocation${DELIM}BioSampleHost${DELIM}BioSampleComment" annot_biosample.txt;
```
Convert from pipe-separated to tab-separated file
```
sed -i "s/|/\t/g" annot_biosample.txt
```
Inspect the annot_biosample.txt file in a spreadsheet view (ex. Excel, Google Sheets)  
Add "REMOVE: Not Yersinia pestis" to the BioSampleComment column to any rows that are confirmed appropriate.  


### Update Database With Annotations
```
nextflow run pipeline.nf \
  --ncbimeta_update ncbimeta.yaml \
  --ncbimeta_annot annot_biosample.txt \
  --max_datasets 3 \
  -resume
```

### Run from established database
```
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets 1 \
  -resume
```

### Tracing
-with-trace
-with-timeline
-with-dag pipeline.pdf
-with-report

### Join Figuring Out
params.ncbimeta_join_first_final = "MasterFirst"
params.ncbimeta_join_first_uniq = "'BioSampleAccession BioSampleAccessionSecondary BioSampleSRAAccession'"
params.ncbimeta_join_first_accessory = "'Assembly SRA'"
params.ncbimeta_join_first_anchor = "BioSample"

// NCBImetaJoin Second Parameters
params.ncbimeta_join_second_final = "Master"
params.ncbimeta_join_second_uniq = "'BioSampleBioProjectAccession'"
params.ncbimeta_join_second_accessory = "'BioProject'"
params.ncbimeta_join_second_anchor = "MasterFirst"

NCBImetaJoin.py \
  --database yersinia_pestis_db.sqlite \
  --anchor BioSample \
  --accessory "Assembly SRA" \
  --final MasterFirst \
  --unique "BioSampleAccession BioSampleAccessionSecondary BioSampleSRAAccession"

NCBImetaJoin.py \
  --database yersinia_pestis_db.sqlite \
  --anchor MasterFirst \
  --accessory BioProject \
  --final Master \

NCBImetaJoin.py \
  --database yersinia_pestis_db.sqlite \
  --anchor BioSample \
  --accessory "Assembly SRA BioProject Nucleotide" \
  --final MasterMinusPubmed \
  --unique "BioSampleAccession BioSampleAccessionSecondary BioSampleSRAAccession BioSampleBioProjectAccession"

NCBImetaJoin.py \
  --database yersinia_pestis_db.sqlite \
  --anchor MasterMinusPubmed \
  --accessory "Pubmed" \
  --final Master \
  --unique "BioSampleAccession BioProjectPublicationID"



### Perl5 Issues
export PERL5LIB=~/miniconda3/envs/phylo-env/lib/site_perl/5.26.2/

### Dev Dependencies
pip install sphinx
pip install sphinx-rtd-theme
pip install m2r
