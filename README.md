# paper-phylogeography
Phylogeography of Yersinia pestis

## Dependencies
**Workflow:** NextFlow  
**Database:** NCBImeta, sqlite3 (CLI)  
**Alignment:** snippy  
**Masking, etc.:** dustmasker, mummer, vcftools   
**Model Selection:** modeltest-ng  
**Statistics:** qualimap, multiqc

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


### Perl5 Issues
```
export PERL5LIB=~/miniconda3/envs/phylo-env/lib/site_perl/5.26.2/
```

### Dev Dependencies
```
pip install sphinx
pip install sphinx-rtd-theme
pip install m2r
```

### SnpEff Build Database
```
mkdir ~/miniconda3/envs/phylo-env/share/snpeff-4.3.1t-3/data/Yersinia_pestis_co92
cp results/reference_genome/GCF_000009065.1_ASM906v1_genomic.gbff ~/miniconda3/envs/phylo-env/share/snpeff-4.3.1t-3/data/Yersinia_pestis_co92/genes.gbk;
# Edit snpEff.config file
Yersinia_pestis_co92.genome : Yersinia_pestis_co92
snpEff build -v -genbank Yersinia_pestis_co92

snpEff -v -csvStats GCA_009669545.1_ASM966954v1_genomic_snippy.snpEff.csv Yersinia_pestis_co92 GCA_009669545.1_ASM966954v1_genomic_snippy.filt.vcf
.filt.vcf
```

### Host selection
```
SELECT BioSampleHost FROM Master WHERE (TRIM(BioSampleHost) > '') AND BioSampleHost IS NOT "missing"
-> 854 Rows
SELECT BioSampleHost FROM Master WHERE (TRIM(BioSampleHost) > '') AND (BioSampleHost IS NOT "missing") AND (BioSampleComment NOT LIKE "%REMOVE%")
-> 734 Rows
SELECT AssemblyFTPGenbank,BioSampleHost,BioSampleCollectionDate,BioSampleGeographicLocation FROM Master WHERE
    (BioSampleComment NOT LIKE "%REMOVE%") AND
	(TRIM(BioSampleHost) > '') AND
	(TRIM(LOWER(BioSampleHost)) IS NOT "missing") AND
	(TRIM(BioSampleCollectionDate) > '') AND
	(TRIM(LOWER(BioSampleCollectionDate)) IS NOT "missing" AND TRIM(LOWER(BioSampleCollectionDate)) IS NOT "not applicable") AND TRIM(LOWER(BioSampleCollectionDate)) IS NOT "unknown" AND TRIM(LOWER(BioSampleCollectionDate)) IS NOT "n/a" ) AND
	(TRIM(BioSampleGeographicLocation) > '') AND
	(TRIM(LOWER(BioSampleCollectionDate)) IS NOT "missing")
->669 Rows
```
