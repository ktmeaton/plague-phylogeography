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

**2. Make a copy of the raw database**

**3. Retain only records relevant to plague or Yersinia pestis.



### EnteroBase
All Yersinia pestis species as of 2018-1207
