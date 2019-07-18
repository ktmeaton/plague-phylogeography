# paper-phylogeography
Phylogeography of Yersinia pestis

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
- Annotate Cui et al. (2013) strains and the Peruvian (2010) strains missing BioProject and/or Strain info
- The misc file fixes problematic strain characters (spaces, parentheses, underscores) as well as adding bioproject info.
```
python NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_cui2013.txt --table BioSample
python NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_peru.txt --table BioSample
python NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_misc.txt --table BioSample

```
- Now keep only the Cui (2013) strains that have BioProject PRJNA47685 and Organization=="BGI". The records that have links to "Beijing Genomics Institute" Organization in this BioProject are duplicates and point to poor BioSample records. They are to be deleted.
- Similarly for the Peruvian project, only the links with Organization=="Institute for Genome Sciences  University of MD School of Medicine" are the good ones with proper links to the SRA. All Peruvian project records with Organization=="University of Maryland School of Medicine, Institute for Genome Sciences" are bad duplicates and should be deleted.
- Remove any BioSample record that does not have BioProject (save for the Black Death project that is "None")

**4. Filter the BioSample table to remove records**
- Mark the BioSample\_id fields with "REMOVE" for undesirable records.
```
python ~/Programs/NCBImeta/src/NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_remove.txt --table BioSample
```
- Records not relevant to plague or Yersinia pestis (ie. Organism/OrganismAlt)
- Transcriptomic sequencing projects ("RNA" in BioSampleTitle or SampleType)
- Laboratory manipulation experiment ("transposon in BioSampleTitle or "vivo" in SampleName or "insertion" in BioSampleTitle)
- Treatment experiment (PRJNA254747 as BioProject, irradiation experiment)
- Prarie dog passage experiments ("prarie dog in BioSampleTitle and PRJNA340278 as BioProject, passage in Infraspecies)
- Duplicate strains (prioritizing original complete genomes, updated version, more complete metadata)

**5. Construct the Master Join table**
```
python ~/Programs/NCBImeta/src/NCBImeta_Join.py --database NCBImeta/output/yersinia_pestis_db.sqlite --anchor BioSample --accessory "BioProject Assembly SRA Nucleotide" --final Master --unique "BioSampleAccession BioSampleAccessionSecondary BioSampleBioProjectAccession"
```



### EnteroBase
All Yersinia pestis species as of 2018-1207
