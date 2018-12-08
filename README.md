# paper2-phylogeography
Phylogeography of Yersinia

## Data Acquisition
Genomic data and metadata acquisition from public repositories.
Desired information falls into 3 categories (sample, project info, data access)  

- Strain
- Alternative Strain Name
- Collection Date
- GeographicFull
- Continent
- Country
- Region
- District
- City
- Latitude
- Longitude
- Elevation
- Environemntal misc (land cover?)
- SourceFull
- Niche
- General or Common Name
- Genus
- Species
- Disease
- Biovar
- Phylogroup
- BioProject Accession
- BioProject Description
- BioSample


Project Information:
- BioProject
- Date entered and modified


Data Access:
- Sequencing Info
- SRA Experiment Run Sample
- Assembly statistics


### NCBI metadata
**1. Retrieve plague metadata from NCBI**

NCBImeta, v0.3.1 commit e53deb1
```
python ~/Programs/NCBImeta/src/NCBImeta.py --flat --config config/paper2-phylogeography_config.py

```

**2. Add columns if they're missing (old program version) vis DB Browser**
```
ALTER TABLE BioSample ADD COLUMN AccessionSecondary TEXT;

ALTER TABLE Assembly ADD COLUMN Comment TEXT;
ALTER TABLE BioSample ADD COLUMN Comment TEXT;
ALTER TABLE BioProject ADD COLUMN Comment TEXT;
ALTER TABLE Nucleotide ADD COLUMN Comment TEXT;
ALTER TABLE SRA ADD COLUMN Comment TEXT;
```

**3. Make a copy of the raw database**
```
cp NCBImeta/output/yersinia_pestis_db.sqlite NCBImeta/output/yersinia_pestis_db_RAW.sqlite
```

**4. Filter the BioSample table to remove:**
- Records not relevant to plague or Yersinia pestis (ie. Organism/OrganismAlt)
- Transcriptomic sequencing projects ("RNA" in BioSampleTitle or SampleType)
- Laboratory manipulation experiment ("transposon in BioSampleTitle or "vivo" in SampleName or "insertion" in BioSampleTitle)
- Treatment experiment (PRJNA254747 as BioProject, irradiation experiment)
- Prarie dog passage experiments ("prarie dog in BioSampleTitle and PRJNA340278 as BioProject, passage in Infraspecies)

**5. Deal with samples that have multiple BioSample accessions or missing BioProject/Strain info**
- Annotate Cui et al. (2013) strains and the Peruvian (2010) strains and misc samples missing BioProject and/or Strain
- The misc script also changes problematic characters (spaces,slashes,brackets)
```
python ~/Programs/NCBImeta/src/NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_cui2013_part2.txt --table BioSample
python ~/Programs/NCBImeta/src/NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_peru_part2.txt --table BioSample
python ~/Programs/NCBImeta/src/NCBImeta_AnnotateReplace.py --database NCBImeta/output/yersinia_pestis_db.sqlite --annotfile NCBImeta/annot/yersinia_pestis_misc.txt --table BioSample

```
- Now keep only the Cui (2013) strains that have BioProject PRJNA47685 and Organization=="BGI". The records that have links to "Beijing Genomics Institute" Organization in this BioProject are broken and to be deleted.
- Similarly for the Peruvian project, only the links with Organization=="University of MD" are the good ones with proper links to the SRA.
- Remove any BioSample record that does not have BioProject (save for the Black Death project that is "None")

**6. Remove duplicate strains**
- Remove duplicate strains
- When considering complete genomes, prefer the original version. Transfer over metadata where needed.
- Prefer assemblies/bioprojects with detailed methodology info.

**7. Construct the Master Join table**



### EnteroBase
All Yersinia pestis species as of 2018-1207
