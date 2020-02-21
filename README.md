# paper-phylogeography
Phylogeography of Yersinia pestis

## Dependencies
NextFlow  
NCBImeta  
sqlite3 (CLI)  
snippy  
dustmasker  
mummer  

## Fresh start

### Build NCBImeta database, test lite run-through of pipeline
```
nextflow run pipeline.nf --ncbimeta_create ncbimeta.yaml --skip_sqlite_import
```

### Remove Wrong Organism Hits
```
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
WHERE (BioSampleOrganism NOT LIKE '%Yersinia pestis%')
```

### Update Database With Annotations
```
nextflow run pipeline.nf \
  --ncbimeta_update ncbimeta.yaml \
  --ncbimeta_annot annot.txt \
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

NCBImetaJoin.py --database yersinia_pestis_db.sqlite --anchor BioSample --accessory "BioProject Assembly SRA Nucleotide" --final Master --unique "BioSampleAccession BioSampleAccessionSecondary BioSampleBioProjectAccession"
