# Plague Phylogeography
A **VERY** in-development work on the phylogeography of *Yersinia pestis*.

## Dependencies
**Workflow:** NextFlow  
**Database:** NCBImeta, sqlite3 (CLI)  
**Alignment:** snippy  
**Masking, etc.:** dustmasker, mummer, vcftools   
**Phylogenetics:** iqtree  
**Statistics:** qualimap, multiqc

### Conda Environment
Create a conda environment with the required dependencies  
```
conda env create -f phylo-env.yaml --name phylo-env
conda activate phylo-env
```

### Dev Dependencies for Docs
```
pip install sphinx sphinx-rtd-theme m2r
```

Everything from here on out is free form notes as I experiment and document.

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
  --max_datasets 2000 \
  -resume
```

### Run from established database
```
nextflow run pipeline.nf \
  --sqlite results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite \
  --max_datasets 1 \
  -resume
```

### Perl5 Issues for HPC Snippy
```
export PERL5LIB=~/miniconda3/envs/phylo-env/lib/site_perl/5.26.2/:$PERL5LIB
```

### SnpEff Build Database
```
mkdir ~/miniconda3/envs/phylo-env/share/snpeff-4.3.1t-3/data/Yersinia_pestis_co92
cp results/reference_genome/GCF_000009065.1_ASM906v1_genomic.gbff ~/miniconda3/envs/phylo-env/share/snpeff-4.3.1t-3/data/Yersinia_pestis_co92/genes.gbk;
# Edit ~/miniconda3/envs/phylo-env/share/snpeff-4.3.1t-3/snpEff.config
Yersinia_pestis_co92.genome : Yersinia_pestis_co92
snpEff build -v -genbank Yersinia_pestis_co92

snpEff -v -csvStats GCA_009669545.1_ASM966954v1_genomic_snippy.snpEff.csv Yersinia_pestis_co92 GCA_009669545.1_ASM966954v1_genomic_snippy.filt.vcf
.filt.vcf
```

### Host selection
Require Assembly,Host,CollectionDate,GeographicLocation
```
SELECT AssemblyFTPGenbank,BioSampleHost,BioSampleCollectionDate,BioSampleGeographicLocation FROM Master WHERE
    (BioSampleComment NOT LIKE "%REMOVE%") AND
	(TRIM(BioSampleHost) > '') AND
	(TRIM(LOWER(BioSampleHost)) IS NOT "missing") AND
	(TRIM(BioSampleCollectionDate) > '') AND
	(TRIM(LOWER(BioSampleCollectionDate)) IS NOT "missing" AND TRIM(LOWER(BioSampleCollectionDate)) IS NOT "not applicable" AND TRIM(LOWER(BioSampleCollectionDate)) IS NOT "unknown" AND TRIM(LOWER(BioSampleCollectionDate)) IS NOT "n/a" ) AND (TRIM(LOWER(BioSampleCollectionDate)) IS NOT "missing") AND
  (TRIM(BioSampleGeographicLocation) > '') AND
  (TRIM(LOWER(BioSampleGeographicLocation)) IS NOT "missing") AND
	(TRIM(AssemblyFTPGenbank) > '')
--> 146 rows
```
Get the "Root" SRA sequences
```
HEADER="Sample_Name\tLibrary_ID\tLane\tSeqType\tOrganism\tStrandedness\tUDG_Treatment\tR1\tR2\tBAM\tGroup\tPopulations\tAge";
echo -e $HEADER > metadata_BronzeAge3.tsv;
sqlite3 results/ncbimeta_db/update/latest/output/database/yersinia_pestis_db.sqlite "SELECT SRASampleName,SRARunAccession, SRALibraryLayout,SRAFileURL From Master WHERE BioSampleAccession IS 'SAMEA3541827'" | \
  awk -F "|" -v Org="Yersinia pestis" -v Strand="double" -v UDG="none" -v BAM="NA" -v Group="NA" -v Pop="NA" -v Age="NA" '{
    Lane = 1; Sample_Name=$1; Library_ID=$2;
    split($2, Library_IDSplit,";");
    split($3, SeqTypeSplit,";");
    split($4, FileURLSplit,";");
    i_new=1
    for (i=1;i<=length(FileURLSplit);i++)
    {
      URL=FileURLSplit[i]
      if (URL ~ /fastq/){FileURLSplitNew[i_new] = URL; i_new++}
    }
    for (i=1;i<=length(Library_IDSplit);i++)
    {
      if (length(SeqTypeSplit) == 1){SeqType = SeqTypeSplit[1]}
      else {SeqType = SeqTypeSplit[i]}
      if (SeqType == "SINGLE"){SeqType = "SE"; R2 = "NA"}
      else if (SeqType == "PAIRED"){SeqType = "PE"; R2 = PLACEHOLDER}
      print Sample_Name "\t" Library_IDSplit[i] "\t" Lane "\t" SeqType "\t" Org "\t" Strand "\t" UDG "\t" FileURLSplitNew[i] "\t" R2 "\t" BAM "\t" Group "\t" Pop "\t" Age;
    }
  }' >> metadata_BronzeAge3.tsv
```


## NCBI Filtering by Attribute AHA!!!
```
"geo_loc_name=Russia: Chechnya"[attr]
```

## NextStrain setup
pip install nextstrain-augur
conda install -c conda-forge --name phylo-env nodejs=10
npm install --global auspice

Activate then deactivate the conda environment to get auspice working!!!

```
mkdir -p results
augur refine \
    --tree ../iqtree/iqtree.core-filter0_bootstrap.treefile \
    --alignment ../snippy_multi/snippy-core.full_CHROM.fasta \
    --vcf-reference ../reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --metadata plagueBasic.tsv \
    --timetree \
    --root residual \
    --coalescent opt \
    --output-tree results/tree.nwk \
    --output-node-data results/branch_lengths.json;
```
```
augur ancestral \
    --tree results/tree.nwk \
    --alignment ../snippy_multi/snippy-core.full.aln \
    --vcf-reference ../reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --inference joint \
    --keep-overhangs \
    --output-node-data results/nt_muts.json
```

```
augur translate \
    --tree results/tree.nwk \
    --ancestral-sequences results/nt_muts.json \
    --vcf-reference ../reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --reference-sequence CHROM.gb \
    --output results/aa_muts.json \
    --alignment-output results/my_alignment_%GENE.fasta \
    --vcf-reference-output results/translations_reference.fasta
```

```
augur traits \
    --tree results/tree.nwk \
    --metadata plagueBasic.tsv \
    --columns region country host \
    --confidence \
    --output results/traits.json
```

Fix the periods (not necessary)
```
sed 's/\.1_/_/g' plagueBasic.tsv > plagueBasic_edit.tsv
sed 's/\.1_/_/g' results/tree.nwk >  results/tree_edit.nwk
for file in `ls results/*json`; do sed 's/\.1_/_/g' $file > ${file%.*}_edit.json; done
```

On export, the output has to be name githubrepo_something
the githubreport and underscore are mandatory! But i think the underscore will screw things up for local deploy...
```
augur export v2 \
    --tree results/tree.nwk \
    --metadata plagueBasic.tsv \
    --node-data results/branch_lengths.json results/traits.json results/aa_muts.json results/nt_muts.json \
    --auspice-config auspice_config.json \
    --output auspice_v2/plague.json
```
Can I make this a smaller export file?? The Default is way too big
```
augur export v2 \
    --tree results/tree.nwk \
    --metadata plagueBasic.tsv \
    --node-data results/branch_lengths.json results/traits.json\
    --auspice-config auspice_config_small.json \
    --lat-longs lat_longs.tsv \
    --output auspice_v2/plague150.json
```
Nicer strain names
```
cp plagueBasic.tsv plagueBasic_edit.tsv;
cp results/tree.nwk results/tree_edit.nwk;
for file in `ls results/*json`; do cp $file ${file%.*}_edit.json; done;
tail -n+2 plagueBasic.tsv | awk -F "\t" '{print $1 "\t" $8}' | while read line;
do
  IN=`echo -e "$line" | cut -f 1`;
  OUT=`echo -e "$line" | cut -f 2`;
  sed -i "s/$IN/$OUT/g" plagueBasic_edit.tsv;
  sed -i "s/$IN/$OUT/g" results/tree_edit.nwk;
  for file in `ls results/*edit.json`; do sed -i "s/$IN/$OUT/g" $file; done;
done
```

```
augur export v2 \
    --tree results/tree_edit.nwk \
    --metadata plagueBasic_edit.tsv \
    --node-data results/branch_lengths_edit.json results/traits_edit.json results/aa_muts_edit.json results/nt_muts_edit.json \
    --auspice-config auspice_config.json \
    --output auspice_v2/plagueEdit.json
```

```
auspice view --datasetDir auspice
```

### Test 150
It looks like the branch support labels from IQTREE and TreeTime are conflicting with each other
Also got a segmentation core dump when running TreeTime large on my laptop.
let's try running IQTREE without branch supports then TreeTime.
(I think branch supports need to go into a json file somehow so they can become --node-data)

```
iqtree \
  -s ../../snippy_multi/snippy-core.full_CHROM.filter0.fasta \
  -m TVM+F+R5 \
  -nt 8 \
  -o Reference \
  -seed 4154541355 \
  -pre iqtree.core-filter0_noBranchSupport \
  2>&1 | tee iqtree.core-filter0_noBranchSupport.output
```

```
mkdir -p results
augur refine \
    --tree iqtree.core-filter0_noBranchSupport.treefile \
    --alignment ../../snippy_multi/snippy-core.full_CHROM.fasta \
    --vcf-reference ../../reference_genome/GCF_000009065.1_ASM906v1_genomic.fna \
    --metadata ../../../nextstrain/plagueBasic.tsv \
    --timetree \
    --root residual \
    --coalescent opt \
    --output-tree results/tree.nwk \
    --output-node-data results/branch_lengths.json;
```

```
augur traits \
    --tree results/tree.nwk \
    --metadata ../../../nextstrain/plagueBasic.tsv \
    --columns region country host \
    --confidence \
    --output results/traits.json
```

```
augur export v2 \
    --tree results/tree.nwk \
    --metadata ../../../nextstrain/plagueBasic.tsv \
    --node-data results/branch_lengths.json results/traits.json \
    --auspice-config ../../../nextstrain/auspice_config_small.json \
    --lat-longs ../../../nextstrain/lat_longs.tsv \
    --output auspice_v2/plague150.json
```
