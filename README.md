# paper-phylogeography
Phylogeography of Yersinia pestis

## Requirements

### Dependencies
NextFlow  
NCBImeta  
sqlite3 (CLI)  
snippy  

## Data Acquisition

### NCBI metadata
**1. Retrieve plague metadata from NCBI**

NCBImeta, v0.3.1 commit 7c10cc8
```
python NCBImeta.py --flat --config NCBImeta/config/paper-phylogeography_config.py

```

**2. Make a copy of the raw database**
```
cp NCBImeta/output/yersinia_pestis_db.sqlite NCBImeta/output/yersinia_pestis_db_RAW.sqlite
```

**3. Deal with samples that have multiple BioSample accessions or missing BioProject/Strain info**
The misc file fixes problematic strain characters (spaces, parentheses, underscores) as well as adding bioproject info.
```
python NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_misc.txt --table BioSample
```
Annotate Cui et al. (2013) strains and the Peruvian (2010) strains missing BioProject and/or Strain info
```
python NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_cui2013.txt --table BioSample
python NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_peru.txt --table BioSample

```

**4. Comment the BioSample table to keep or remove records**
Mark the BioSampleComment fields with "REMOVE" for undesirable records.

```
python ~/Programs/NCBImeta/src/NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_REMOVE.txt --table BioSample
```
- Records not relevant to plague or Yersinia pestis (ie. Organism/OrganismAlt)
- Transcriptomic sequencing projects (ex. "RNA" in BioSampleTitle or SampleType)
- Laboratory manipulation experiment (ex. "transposon in BioSampleTitle or "vivo" in SampleName or "insertion" in BioSampleTitle)
- Treatment experiment (PRJNA254747 as BioProject, irradiation experiment)
- Prarie dog passage experiments ("prarie dog in BioSampleTitle and PRJNA340278 as BioProject, passage in Infraspecies)
- [?] Duplicate strains (prioritizing original complete genomes, updated version, more complete metadata)

Mark the BioSampleComment fields with "KEEP" for desirable records.  
```
python ~/Programs/NCBImeta/src/NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_KEEP.txt --table BioSample
```

Any record that now has an empty comment field has been added to the database after SEPT-2019.  



**5. Construct the Master Join table**
```
python ~/Programs/NCBImeta/src/NCBImeta_Join.py --database NCBImeta/output/yersinia_pestis_db.sqlite --anchor BioSample --accessory "BioProject Assembly SRA Nucleotide" --final Master --unique "BioSampleAccession BioSampleAccessionSecondary BioSampleBioProjectAccession"
```



### EnteroBase
All Yersinia pestis species as of 2018-1207
