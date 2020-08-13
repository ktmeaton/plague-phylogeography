
Local Reads Prep
----------------

Prepare custom read data as eager tsv.

**shell**::

	biosampleColumn=1
	inTSV=!{params.eager_tsv}
	outTSV=`basename ${inTSV%.*}.txt`
	tail -n+2 !{custom_tsv_eager} | cut -f ${biosampleColumn} | sort | uniq > ${outTSV}

Local Assembly Prep
-------------------

Prepare custom assembly data as file list.


Ncbimeta Db Create
------------------

Run NCBImeta queries to generate db from scratch.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_ncbimeta_yaml                         yaml                                     NCBImeta config file.                    
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_ncbimeta_sqlite_update                sqlite                                   NCBImeta SQLite database for process :ref:`ncbimeta_db_update<Ncbimeta_Db_Update>`
ch_ncbimeta_yaml_update                  yaml                                     NCBImeta config file for process :ref:`ncbimeta_db_update<Ncbimeta_Db_Update>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
ncbimeta_sqlite_db                       sqlite                                   NCBImeta SQLite database.                
ncbimeta_yaml                            yaml                                     NCBImeta config file.                    
*.log                                    text                                     Text logs of NCBImeta database creation. 
======================================== ======================================== ========================================

**script**::

	NCBImeta.py --config ${ncbimeta_yaml}

Ncbimeta Db Update
------------------

Run NCBImeta queries to update, annotate, and join a previously created database.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_ncbimeta_yaml_update                  yaml                                     NCBImeta config file from process :ref:`ncbimeta_db_create<Ncbimeta_Db_Create>`
ch_ncbimeta_sqlite_update                sqlite                                   NCBImeta SQLite database from process :ref:`ncbimeta_db_create<Ncbimeta_Db_Create>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_ncbimeta_sqlite_import                sqlite                                   NCBImeta SQLite database for process :ref:`sqlite_import<Sqlite_Import>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
ncbimeta_yaml                            yaml                                     NCBImeta config file.                    
*.log                                    text                                     Text logs of NCBImeta database update.   
*.txt                                    text                                     Text export of NCBImeta database.        
======================================== ======================================== ========================================

**script**::

	# Make directories to mirror NCBImeta expected structure
	mkdir ${params.ncbimeta_output_dir};
	mkdir ${params.ncbimeta_output_dir}/database;
	mkdir ${params.ncbimeta_output_dir}/log;
	# Copy over input files
	cp ${ncbimeta_sqlite} ${params.ncbimeta_output_dir}/database;
	cp ${outdir}/ncbimeta_db/update/latest/${params.ncbimeta_output_dir}/log/* ${params.ncbimeta_output_dir}/log;
	# Execute NCBImeta
	NCBImeta.py --config ${ncbimeta_yaml}
	# If annotation file supplied, run the annotation script
	if [[ ${params.ncbimeta_annot} != "false" ]]; then
	ANNOT_FILE=`basename ${params.ncbimeta_annot}`
	mv ${workDir}/dummy_annot.txt `pwd`/\$ANNOT_FILE;
	NCBImetaAnnotateReplace.py --table ${params.ncbimeta_annot_table} --annot ${params.ncbimeta_annot} --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db}
	fi
	# Drop old or outdated join tables
	sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS MasterFirst"
	sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS MasterSecond"
	sqlite3 ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} "DROP TABLE IF EXISTS Master"
	# Join Tables
	NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_first_anchor} --accessory ${params.ncbimeta_join_first_accessory} --final ${params.ncbimeta_join_first_final} --unique ${params.ncbimeta_join_first_uniq}
	NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_second_anchor} --accessory ${params.ncbimeta_join_second_accessory} --final ${params.ncbimeta_join_second_final} --unique ${params.ncbimeta_join_second_uniq}
	NCBImetaJoin.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --anchor ${params.ncbimeta_join_third_anchor} --accessory ${params.ncbimeta_join_third_accessory} --final ${params.ncbimeta_join_third_final} --unique ${params.ncbimeta_join_third_uniq}
	# Export Tables
	NCBImetaExport.py --database ${params.ncbimeta_output_dir}/database/${params.ncbimeta_sqlite_db} --outputdir ${params.ncbimeta_output_dir}/database/

Sqlite Import
-------------

Import assembly FTP url from database, retrieve file names for web get, prepare TSV input of SRA metadata for EAGER pipeline.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_sqlite                                sqlite                                   NCBImeta SQLite database from process :ref:`ncbimeta_db_update or params.sqlite<Ncbimeta_Db_Update Or Params.Sqlite>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_assembly_download_ftp                 text                                     FTP url for process :ref:`assembly_download<Assembly_Download>`
ch_sra_tsv_eager                         tsv                                      TSV metadata input for process :ref:`eager<Eager>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
file_assembly_download_ftp               text                                     List of FTP urls for genomic assembly download.
eager_tsv                                tsv                                      TSV metadata input for EAGER pipeline.   
======================================== ======================================== ========================================

**shell**::

	# Select the Genbank Assemblies
	if [[ !{params.sqlite_select_command_asm} != "false"  ]]; then
	sqlite3 !{sqlite} !{params.sqlite_select_command_asm} | \
	grep . | \
	head -n !{params.max_datasets_assembly} | \
	sed -E -e 's/ |;/\\n/g' | \
	while read line;
	do
	if [[ ! -z ${line} ]]; then
	asm_ftp=`echo ${line} | \
	awk -F "/" -v suffix=!{params.genbank_assembly_gz_suffix} '{print $0 FS $NF suffix}'`;
	echo ${asm_ftp} >> !{params.file_assembly_download_ftp}
	fi;
	done;
	fi;
	
	# Extract SRA Metadata for EAGER tsv
	if [[ !{params.sqlite_select_command_sra} != "false"  ]]; then
	!{params.scriptdir}/sqlite_EAGER_tsv.py \
	--database !{sqlite} \
	--query !{params.sqlite_select_command_sra} \
	--organism !{params.eager_organism} \
	--max-datasets !{params.max_datasets_sra} \
	--output metadata_sra_eager.tsv \
	--fastq-dir !{outdir}/sra_download/
	biosampleColumn=1
	accessionColumn=2
	tail -n+2 metadata_sra_eager.tsv | cut -f $biosampleColumn | sort | uniq > metadata_sra_biosample.tsv
	fi;

Assembly Download
-----------------

Download genomic assembly fasta using FTP urls.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_assembly_fna_gz_local                 fasta.gz                                 The genomic assembly accessed by url via FTP.
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_assembly_fna_snippy_pairwise          fasta                                    The genomic assembly for process :ref:`snippy_pairwise<Snippy_Pairwise>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
genbank_assembly_fna_suffix              fasta                                    The locally downloaded genomic assembly. 
======================================== ======================================== ========================================

**script**::

	# Use -f otherwise error due to too many levels of symbolic links
	gunzip -f ${assembly_fna_gz}

Sra Download
------------

Download sequence data from the SRA database.

**script**::

	# Change the download sra location and timeout settings
	mkdir -p ~/.ncbi/
	# Default sra cache path
	sra_fastq_dump_path=${sra_fastq_dump_path}
	
	# Create SRA config file if it doesn't exist
	if [[ ! -f $HOME/.ncbi/user-settings.mkfg ]]; then
	echo '/repository/user/main/public/root = "\${sra_fastq_dump_path}"' > $HOME/.ncbi/user-settings.mkfg
	fi
	
	# Set cache enabled if not set
	if [[ -z `grep "/cache-enabled" $HOME/.ncbi/user-settings.mkfg` ]]; then
	echo '/cache-enabled = "true"' >> $HOME/.ncbi/user-settings.mkfg
	fi;
	
	# Set the cache path
	if [[ -z `grep "/repository/user/main/public/root" $HOME/.ncbi/user-settings.mkfg` ]]; then\
	# Set SRA Cache Path
	echo '/repository/user/main/public/root = "\${sra_fastq_dump_path}"' >> $HOME/.ncbi/user-settings.mkfg
	else
	# Retrieve SRA Cache Path
	sra_fastq_dump_path=`grep "/repository/user/main/public/root" $HOME/.ncbi/user-settings.mkfg | \
	cut -d " " -f 3 | \
	sed 's/"//g'`
	fi;
	
	# Set the timeout
	if [[ -z `grep "/http/timeout/read" $HOME/.ncbi/user-settings.mkfg` ]]; then
	echo '/http/timeout/read = "10000"' >> $HOME/.ncbi/user-settings.mkfg
	fi;
	
	echo "SRA Cache:" \${sra_fastq_dump_path}
	echo "NCBI settings:" `cat $HOME/.ncbi/user-settings.mkfg`
	
	# Create organization directories
	mkdir -p ${sra_biosample_val}
	mkdir -p ${sra_biosample_val}/single;
	mkdir -p ${sra_biosample_val}/paired;
	
	# Retrieve sra accessions for the biosample
	accessionCol=2
	sraAccList=`grep -w ${sra_biosample_val} ${tsv_eager} | cut -f \$accessionCol`;
	for sraAcc in \$sraAccList;
	do
	validate='false'
	# Keep trying to download until valid file is acquired
	while [ \$validate == 'false' ]
	do
	# Download fastq files from the SRA
	fastq-dump \
	--outdir ${sra_biosample_val}/ \
	--skip-technical \
	--gzip \
	--split-files \$sraAcc;
	# Validate sra file
	ls -l \${sra_fastq_dump_path}/sra/\${sraAcc}.sra*
	validate_str=`vdb-validate \${sra_fastq_dump_path}/sra/\${sraAcc}.sra* 2>&1`
	echo \${validate_str}
	if [[ \${validate_str} != *"corrupt"* ]]; then
	validate='true'
	else
	echo "Removing \${sraAcc} from the SRA cache."
	rm \${sra_fastq_dump_path}/sra/\${sraAcc}.sra*
	fi
	done
	
	# If a paired-end or single-end file was downloaded
	if [ -f ${sra_biosample_val}/\${sraAcc}_1.fastq.gz ] &&
	[ -f ${sra_biosample_val}/\${sraAcc}_2.fastq.gz ]; then
	mv ${sra_biosample_val}/\${sraAcc}*.fastq.gz ${sra_biosample_val}/paired/;
	else
	mv ${sra_biosample_val}/\${sraAcc}*.fastq.gz ${sra_biosample_val}/single/;
	fi
	done

Reference Download
------------------

Download the reference genome of interest from the FTP site.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
reference_genome_fna_ftp                 fasta.gz                                 The reference genome fasta accessed by url via FTP.
reference_genome_gb_ftp                  fasta.gz                                 The reference genome gbff accessed by url via FTP.
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_reference_detect_repeats              fasta                                    The reference genome for process :ref:`detect_repeats<Detect_Repeats>`
ch_reference_genome_detect_low_complexity fasta                                    The reference genome for process :ref:`detect_low_complexity<Detect_Low_Complexity>`
ch_reference_gb_snippy_pairwise          gbff                                     The reference genome for process :ref:`snippy_pairwise<Snippy_Pairwise>`
ch_reference_gb_snippy_multi             gbff                                     The reference genome for process :ref:`snippy_multi<Snippy_Multi>`
ch_reference_genome_snpeff_build_db      gbff                                     The reference genome for process :ref:`snpeff_build_db<Snpeff_Build_Db>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
reference_genome_fna_local               fasta                                    The locally downloaded reference fasta.  
reference_genome_gb_local                gbff                                     The locally downloaded reference annotations.
======================================== ======================================== ========================================

**script**::

	gunzip -f ${reference_genome_fna_local}
	gunzip -f ${reference_genome_gb_local}
	gunzip -f ${reference_genome_gff_local}
	# Edit the fasta headers to match the gb loci (for snippy)
	GB_LOCI=(`grep LOCUS ${reference_genome_gb_local.baseName} | sed 's/ \\+/ /g' | cut -d " " -f 2`);
	FNA_LOCI=(`grep ">" ${reference_genome_fna_local.baseName} | cut -d " " -f 1 | cut -d ">" -f 2`);
	i=0;
	while [ \$i -lt \${#GB_LOCI[*]} ];
	do
	sed -i "s/\${FNA_LOCI[\$i]}/\${GB_LOCI[\$i]}/g" ${reference_genome_fna_local.baseName};
	i=\$(( \$i + 1));
	done
	# Extract chromosome sequence
	CHROM=NC_003143
	fnaName=${reference_genome_fna_local.baseName}
	fnaNameCHROM=\${fnaName%.*}_CHROM.fna
	samtools faidx ${reference_genome_fna_local.baseName};
	samtools faidx ${reference_genome_fna_local.baseName} \${CHROM} \
	> \$fnaNameCHROM
	

Snpeff Build Db
---------------

Build a SnpEff database for the reference genome annotations.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
reference_genome_gb                      gbff                                     The reference genome gbff from process :ref:`reference_download<Reference_Download>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snpeff_config_snippy_pairwise         text                                     Edited SnpEff configuration file for process :ref:`snippy_pairwise<Snippy_Pairwise>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
snpEff.config                            text                                     Edited SnpEff configuration file.        
snpEffectPredictor.bin                   gzip text                                SnpEff database.                         
======================================== ======================================== ========================================

**script**::

	# Locate SnpEff directories in miniconda path
	ref=${reference_genome_gb.baseName}
	snpeffDir=\${CONDA_PREFIX}/share/snpeff*
	snpeffData=\$snpeffDir/data;
	
	# Make a SnpEff database dir
	mkdir -p data/
	mkdir -p data/\$ref/
	
	# Move over the reference genbank annotations and rename
	cp ${reference_genome_gb} data/\$ref/genes.gbk;
	
	# Copy over snpEff.config
	cp \$snpeffDir/snpEff.config .
	
	# Add the new annotation entry to the snpeff config file
	configLine="${reference_genome_gb.baseName}.genome : ${reference_genome_gb.baseName}"
	
	# Search for the genome entry in the snpEff config file
	if [[ -z `grep "\$configLine" snpEff.config` ]]; then
	echo "\$configLine" >> snpEff.config;
	fi;
	
	# Build the snpEff databse
	snpEff build -dataDir ./data/ -v -genbank ${reference_genome_gb.baseName}

Reference Detect Repeats
------------------------

Detect in-exact repeats in reference genome with mummer and convert the identified regions file to bed format.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_reference_genome_detect_repeats       fasta                                    The reference genome fasta from the process :ref:`reference_download<Reference_Download>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_bed_ref_detect_repeats                bed                                      A bed file containing regions of in-exact repeats for process :ref:`snippy_merge_mask_bed<Snippy_Merge_Mask_Bed>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
reference_genome_fna.inexact.coords      coords                                   Alignment coordinate file generated by mummer.
reference_genome_fna.inexact.repeats     coords                                   Filtered file for sequence similarity and self-alignments
reference_genome_fna.inexact.repeats.bed bed                                      Bed file created from filtered coordinates and adjusted for 0-base system.
======================================== ======================================== ========================================

**script**::

	PREFIX=${reference_genome_fna.baseName}
	# Align reference to itself to find inexact repeats
	nucmer --maxmatch --nosimplify --prefix=\${PREFIX}.inexact ${reference_genome_fna} ${reference_genome_fna}
	# Convert the delta file to a simplified, tab-delimited coordinate file
	show-coords -r -c -l -T \${PREFIX}.inexact.delta | tail -n+5 > \${PREFIX}.inexact.coords
	# Remove all "repeats" that are simply each reference aligned to itself
	# also retain only repeats with more than 90% sequence similarity.
	awk -F "\t" '{if (\$1 == \$3 && \$2 == \$4 && \$12 == \$13)
	{next;}
	else if (\$7 > 90)
	{print \$0}}' \${PREFIX}.inexact.coords > \${PREFIX}.inexact.repeats
	# Also exact and tandem repeats??
	# Convert to bed file format, changing to 0-base position coordinates
	awk -F "\t" '{print \$12 "\t" \$1-1 "\t" \$2-1;
	if (\$3 > \$4){tmp=\$4; \$4=\$3; \$3=tmp;}
	print \$13 "\t" \$3-1 "\t" \$4-1;}' \${PREFIX}.inexact.repeats | \
	sort -k1,1 -k2,2n | \
	bedtools merge > \${PREFIX}.inexact.repeats.bed

Reference Detect Low Complexity
-------------------------------

Detect low complexity regions with dustmasker and convert the identified regions file to bed format.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_reference_genome_low_complexity       fasta                                    The reference genome fasta from the process :ref:`reference_download<Reference_Download>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_bed_ref_low_complex                   bed                                      A bed file containing regions of low-complexity regions for process :ref:`snippy_merge_mask_bed<Snippy_Merge_Mask_Bed>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
reference_genome_fna.dustmasker.intervals intervals                                Interval file containing regions of low-complexity.
reference_genome_fna.dustmasker.bed      bed                                      Bed file created from intervals and adjusted for 0-base system.
======================================== ======================================== ========================================

**script**::

	dustmasker -in ${reference_genome_fna} -outfmt interval > ${reference_genome_fna.baseName}.dustmasker.intervals
	${params.scriptdir}/intervals2bed.sh ${reference_genome_fna.baseName}.dustmasker.intervals ${reference_genome_fna.baseName}.dustmasker.bed

Outgroup Download
-----------------

Download the outgroup assemblies.

**script**::

	gunzip -f ${outgroup_fna_local}
	# Store the file basename/prefix for iqtree outgroup param
	filename=${outgroup_fna_local}
	fna="\${filename%.*}"
	prefix="\${fna%.*}"

Eager
-----

Run the nf-core/eager pipeline on SRA samples.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_reference_genome_eager                fna                                      The reference genome fasta from process :ref:`reference_genome_download<Reference_Genome_Download>`
ch_sra_fastq_eager                       fastq                                    The sra fastq sequences from process :ref:`sra_download<Sra_Download>`
ch_tsv_eager                             tsv                                      The sra metadata tsv from process :ref:`sqlite_import<Sqlite_Import>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_sra_bam_snippy_pairwise               fastq                                    The deduplicated aligned bam for process :ref:`snippy_pairwise<Snippy_Pairwise>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
damageprofiler/*                         misc                                     aDNA damage visualization and statistics.
deduplication/*                          misc                                     Deduplicated aligned bam and statistics. 
pipeline_info/*                          misc                                     Pipeline information.                    
preseq/*                                 misc                                     Preseq complexity statistics.            
qualimap/*                               misc                                     Genome coverage and depth visualization and statistics.
MultiQC/*                                misc                                     Multi software visualizations and statistics.
SoftwareVersions/*                       misc                                     Version of all software used in nf-core eager.
======================================== ======================================== ========================================

**shell**::

	# Create biosample specific tsv input for eager
	head -n 1 !{eager_tsv} > metadata_!{biosample_val}.tsv
	grep -w !{biosample_val} !{eager_tsv} >> metadata_!{biosample_val}.tsv
	
	# The set command is to deal with PS1 errors
	set +eu
	# Enable conda activate support in this bash subshell
	CONDA_BASE=$(conda info --base) ;
	source ${CONDA_BASE}/etc/profile.d/conda.sh
	
	# Activate the eager environment
	conda activate nf-core-eager-2.2.0dev
	
	# Run the eager command
	task_mem_reformat=`echo !{task.memory} | sed 's/ /./g'`
	nextflow -C ~/.nextflow/assets/nf-core/eager/nextflow.config \
	run nf-core/eager \
	-r !{params.eager_rev} \
	--input metadata_!{biosample_val}.tsv \
	--outdir . \
	--fasta !{reference_genome_fna} \
	--clip_readlength !{params.eager_clip_readlength} \
	--preserve5p \
	--mergedonly \
	--mapper bwaaln \
	--bwaalnn !{params.eager_bwaalnn} \
	--bwaalnl !{params.eager_bwaalnl} \
	--run_bam_filtering \
	--bam_mapping_quality_threshold !{params.snippy_map_qual} \
	--bam_discard_unmapped \
	--bam_unmapped_type discard \
	--max_memory ${task_mem_reformat} \
	--max_cpus !{task.cpus} \
	--max_time !{task.time}
	
	# Deactivate the eager env
	conda deactivate
	set +eu
	
	# Rename deduplication bam for snippy pairwise RG
	dir="final_bams"
	mkdir -p $dir;
	if [[ -d merged_bams/ ]]; then
	mergedBam=`ls merged_bams/*/*.bam`;
	else
	mergedBam=`ls deduplication/*/*.bam`;
	fi
	for file in `ls ${mergedBam}`;
	do
	outfile=$dir/!{biosample_val}.bam;
	samtools addreplacerg -r ID:!{biosample_val} -r SM:!{biosample_val} -o $outfile $file
	done
	
	# Move pipeline trace and multiqc into named sample folder
	mkdir -p pipeline_info/!{biosample_val}/
	mv pipeline_info/*txt pipeline_info/*html pipeline_info/*svg pipeline_info/!{biosample_val}/
	mkdir -p MultiQC/!{biosample_val}/
	mv MultiQC/multiqc_data/ MultiQC/multiqc_report.html MultiQC/!{biosample_val}/

Snippy Pairwise
---------------

Pairwise align contigs to reference genome with snippy.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_assembly_fna_snippy_pairwise          fasta                                    The genomic assembly from process :ref:`assembly_download<Assembly_Download>`
ch_reference_gb_snippy_pairwise          gbff                                     The reference annotations from process :ref:`reference_download<Reference_Download>`
ch_snpeff_config_snippy_pairwise         text                                     Edited SnpEff configuration file from process :ref:`snpeff_build_db<Snpeff_Build_Db>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_snps_variant_summary           text                                     Table of summarized SNP counts for process :ref:`variant_summary<Variant_Summary>`
ch_snippy_subs_vcf_detect_density        vcf                                      Substitutions for process :ref:`pairwise_detect_snp_high_density<Pairwise_Detect_Snp_High_Density>`
ch_snippy_bam_pairwise_qualimap          bam                                      Pairwise alignment file for process :ref:`qualimap_snippy_pairwise<Qualimap_Snippy_Pairwise>`
ch_snippy_csv_snpEff_multiqc             csv                                      Variant summary statistics for process :ref:`multiqc<Multiqc>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
assembly_fna_snippy.summary.txt          text                                     Table of summarized SNP counts.          
assembly_fna_snippy.subs.vcf             vcf                                      Substitutions.                           
assembly_fna_snippy.csv                  csv                                      SnpEff annotation and summary report.    
assembly_fna_snippy.bam                  bam                                      Snippy bam alignment file.               
assembly_fna_snippy.*                    misc                                     All default snippy pipeline output.      
======================================== ======================================== ========================================

**script**::

	if [[ "${fna_bam.extension}" == "fna" ]]; then
	snippy \
	--prefix ${fna_bam.baseName}_snippy \
	--cpus ${task.cpus} \
	--reference ${reference_genome_gb} \
	--outdir output${params.snippy_ctg_depth}X/${fna_bam.baseName} \
	--ctgs ${fna_bam} \
	--mapqual ${params.snippy_map_qual} \
	--mincov ${params.snippy_ctg_depth} \
	--minfrac ${params.snippy_min_frac} \
	--basequal ${params.snippy_base_qual} \
	--report;
	elif  [[ "${fna_bam.extension}" == "bam" ]]; then
	snippy \
	--prefix ${fna_bam.baseName}_snippy \
	--cpus ${task.cpus} \
	--reference ${reference_genome_gb} \
	--outdir output${params.snippy_ctg_depth}X/${fna_bam.baseName} \
	--bam ${fna_bam} \
	--mapqual ${params.snippy_map_qual} \
	--mincov ${params.snippy_bam_depth} \
	--minfrac ${params.snippy_min_frac} \
	--basequal ${params.snippy_base_qual} \
	--report;
	fi;
	
	# Save Output Dir for snippy_multi channel
	snippyDir=`pwd`"/output${params.snippy_ctg_depth}X/${fna_bam.baseName}/"
	
	snippy_snps_in=output${params.snippy_ctg_depth}X/${fna_bam.baseName}/${fna_bam.baseName}_snippy.txt
	snippy_snps_txt=output${params.snippy_ctg_depth}X/${fna_bam.baseName}/${fna_bam.baseName}_snippy.summary.txt
	
	COMPLEX=`awk 'BEGIN{count=0}{if (\$1 == "Variant-COMPLEX"){count=\$2}}END{print count}' \$snippy_snps_in;`
	DEL=`awk 'BEGIN{count=0}{if (\$1 == "Variant-DEL"){count=\$2}}END{print count}' \$snippy_snps_in;`
	INS=`awk 'BEGIN{count=0}{if (\$1 == "Variant-INS"){count=\$2}}END{print count}' \$snippy_snps_in;`
	MNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-MNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
	SNP=`awk 'BEGIN{count=0}{if (\$1 == "Variant-SNP"){count=\$2}}END{print count}' \$snippy_snps_in;`
	TOTAL=`awk 'BEGIN{count=0}{if (\$1 == "VariantTotal"){count=\$2}}END{print count}' \$snippy_snps_in;`
	echo -e output${params.snippy_ctg_depth}X/${fna_bam.baseName}"\\t"\$COMPLEX"\\t"\$DEL"\\t"\$INS"\\t"\$MNP"\\t"\$SNP"\\t"\$TOTAL >> \$snippy_snps_txt
	
	snippy_snps_filt=output${params.snippy_ctg_depth}X/${fna_bam.baseName}/${fna_bam.baseName}_snippy.filt.vcf
	snippy_snps_csv=output${params.snippy_ctg_depth}X/${fna_bam.baseName}/${fna_bam.baseName}_snippy.csv
	snippy_snps_rename=output${params.snippy_ctg_depth}X/${fna_bam.baseName}/${fna_bam.baseName}_snippy.rename.csv
	
	# SnpEff csv Stats
	mv \$snippy_snps_csv \$snippy_snps_rename
	snpEff -c ${snpeff_config} \
	-dataDir ${outdir}/reference_genome/data/ \
	-csvStats \$snippy_snps_csv \
	-quiet \
	${reference_genome_gb.baseName} \
	\$snippy_snps_filt

Snippy Variant Summary Collect
------------------------------

Concatenate variant summary tables for all samples.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_snps_variant_summary           text                                     Table of single-sample summarized SNP counts from process :ref:`snippy_pairwise<Snippy_Pairwise>`
ch_snippy_variant_summary_multi_collect  text                                     Table of multi-sample summarized SNP counts.
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_variant_summary_multiqc        text                                     Table of multi-sample summarized SNP counts for process :ref:`multiqc<Multiqc>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
snippy_variant_summary.txt               text                                     Table of multi-sample summarized SNP counts.
======================================== ======================================== ========================================


Snippy Detect Snp High Density
------------------------------

Detect regions of high SNP density.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_subs_vcf_detect_density        vcf                                      Substitutions from process :ref:`snippy_pairwise<Snippy_Pairwise>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_subs_bed_merge_density         bed                                      High-density SNP regions for process :ref:`snippy_merge_snp_high_density<Snippy_Merge_Snp_High_Density>`
======================================== ======================================== ========================================

**script**::

	vcftools --vcf ${snippy_subs_vcf} --SNPdensity ${params.snippy_snp_density_window} --out ${snippy_subs_vcf.baseName}.tmp
	tail -n+2 ${snippy_subs_vcf.baseName}.tmp.snpden | awk -F "\\t" '{if (\$3 > 1){print \$1 "\\t" \$2-10-1 "\\t" \$2}}' > ${snippy_subs_vcf.baseName}.snpden

Snippy Sort Snp High Density
----------------------------

Sort and merge regions of high SNP density.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_subs_bed_sort_density          bed                                      High density SNP regions collected after process :ref:`snippy_detect_snp_high_density<Snippy_Detect_Snp_High_Density>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_subs_bed_density_multi         bed                                      Sorted and merged high density SNP regions for process :ref:`snippy_multi<Snippy_Multi>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
snippy_variant_density                   bed                                      Sorted and merged high density SNP regions.
======================================== ======================================== ========================================

**script**::

	sort -k1,1 -k2,2n ${snippy_subs_bed} | bedtools merge > ${params.snippy_variant_density}.txt

Snippy Merge Mask Bed
---------------------

Combine, merge, and sort all BED file regions for masking the multiple alignment.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_bed_ref_detect_repeats                bed                                      A bed file containing regions of in-exact repeats from process :ref:`reference_detect_repeats<Reference_Detect_Repeats>`
ch_bed_ref_low_complex                   bed                                      A bed file containing regions of low-complexity regions from process :ref:`reference_detect_low_complexity<Reference_Detect_Low_Complexity>`
ch_snippy_subs_bed_density_multi         bed                                      Sorted and merged high density SNP regions from process :ref:`snippy_sort_snp_high_density<Snippy_Sort_Snp_High_Density>`
ch_bed_mask_master_merge                 bed                                      Combined BED files of repeats, low-complexity and 
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_bed_mask_snippy_multi                 bed                                      Master masking BED file for process :ref:`snippy_multi<Snippy_Multi>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
master.bed                               bed                                      Master masking BED file.                 
======================================== ======================================== ========================================

**script**::

	cat ${bed_mask} | sort -k1,1 -k2,2n | bedtools merge > master.bed

Snippy Multi
------------

Perform a multiple genome alignment with snippy-core.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_reference_gb_snippy_multi             gbff                                     The reference genome from process :ref:`reference_download<Reference_Download>`
ch_bed_mask_snippy_multi                 bed                                      Master masking BED file from process :ref:`snippy_merge_mask_bed<Snippy_Merge_Mask_Bed>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_core_aln_filter                fasta                                    Multi fasta of aligned core SNPs for process :ref:`snippy_multi_filter<Snippy_Multi_Filter>`
ch_snippy_core_full_aln_filter           fasta                                    Multi fasta of aligned core genome for process :ref:`snippy_multi_filter<Snippy_Multi_Filter>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
\*                                       misc                                     All default output from snippy-core.     
======================================== ======================================== ========================================

**script**::

	# Store a list of all the Snippy output directories in a file
	allDir=`for path in ${snippy_outdir_path};
	do
	echo \$path | sed 's/\\[\\|,\\|\\]//g' ;
	done | tr '\n' ' ' `;
	
	# Perform multiple genome alignment (with custom filtering)
	snippy-core \
	--ref ${reference_genome_gb} \
	--prefix snippy-core \
	--mask ${bed_mask} \
	--mask-char ${params.snippy_mask_char} \
	\$allDir 2>&1 | tee snippy-core.log

Snippy Multi Filter
-------------------

Filter the multiple alignment for X% missing data and split by locus.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_core_full_aln_filter           fasta                                    Multi fasta of aligned core genome ffrom process :ref:`snippy_multi<Snippy_Multi>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_core_filter_iqtree             fasta                                    Multi fasta of filtered core genome sites for process :ref:`iqtree<Iqtree>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
snippy_core_full_aln.filter\*.fasta      fasta                                    Multi fasta of filtered chromosome genome sites.
*.fasta                                  fasta                                    All loci extracted fasta files.          
*.bed                                    bed                                      All loci bed coordinate files for extraction.
======================================== ======================================== ========================================

**script**::

	# Split by LOCUS (generates snippy-core_%REPLICON.fasta)
	${params.scriptdir}/fasta_split_locus.sh ${snippy_core_full_aln}
	# Filter full CHROMOSOME alignment (No Missing Data)
	snp-sites -m -c -b -o ${snippy_core_full_aln.baseName}_CHROM.filter0.fasta ${snippy_core_full_aln.baseName}_CHROM.fasta;
	# Optional: Filter full alignment to remove less missing data
	if [[ ${params.snippy_multi_missing_data_text} > 0 ]]; then
	${params.scriptdir}/fasta_unwrap.sh ${snippy_core_full_aln.baseName}_CHROM.fasta > ${snippy_core_full_aln.baseName}_CHROM.unwrap.fasta;
	${params.scriptdir}/fasta_filterGapsNs.sh \
	${snippy_core_full_aln.baseName}_CHROM.unwrap.fasta \
	${params.snippy_multi_missing_data} \
	${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.backbone > \
	${snippy_core_full_aln.baseName}_CHROM.filter${params.snippy_multi_missing_data_text}.fasta;
	fi;

Iqtree
------

Maximum likelihood tree search and model selection, iqtree phylogeny.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_core_filter_iqtree             fasta                                    Multi fasta of filtered core genome sites from process :ref:`snippy_multi_filter<Snippy_Multi_Filter>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_iqtree_treefile_augur_refine          newick                                   Newick treefile phylogeny with branch supports for process :ref:`augur_refine<Augur_Refine>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
iqtree.core-filter*_bootstrap.treefile   newick                                   Newick treefile phylogeny with branch supports.
!*treefile                               misc                                     All default output of iqtree other than the treefile.
======================================== ======================================== ========================================

**script**::

	# Setup the outgroup
	if [[ ${params.skip_outgroup_download} == "false"  ]]; then
	OUTGROUP="${outgroup_file}";
	# Strip brackets and spaces from list
	OUTGROUP=`echo "\$OUTGROUP" | sed 's/\\[\\| \\|\\]//g'`;
	else
	OUTGROUP=${params.iqtree_outgroup}
	fi
	
	# Setup the model or model testing
	if [[ ${params.iqtree_model} == "false"  ]]; then
	MODEL="MFP"
	else
	MODEL="${params.iqtree_model}"
	fi
	
	# Setup the branch support param
	if [[ ${params.iqtree_branch_support} == "true"  ]]; then
	BRANCH_SUPPORT="--ufboot ${params.iqtree_ufboot} --alrt ${params.iqtree_ufboot}";
	SUFFIX="_bootstrap";
	else
	BRANCH_SUPPORT="";
	SUFFIX="";
	fi
	
	# A thorough tree search for model selection can be done with -m MF -mtree
	iqtree \
	-s ${snippy_core_filter_aln} \
	-m \$MODEL \
	--threads-max ${task.cpus} \
	-nt AUTO \
	-o \$OUTGROUP \
	-seed \$RANDOM \
	\${BRANCH_SUPPORT} \
	--runs ${params.iqtree_runs} \
	-pre iqtree.core-filter${params.snippy_multi_missing_data_text}\${SUFFIX} \
	2>&1 | tee iqtree.core-filter${params.snippy_multi_missing_data_text}\${SUFFIX}.output

Qualimap Snippy Pairwise
------------------------

Run QualiMap on the output bam of snippy pairwise.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_bam_pairwise_qualimap          bam                                      Pairwise alignment file from process :ref:`snippy_pairwise<Snippy_Pairwise>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Output:                                  Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_pairwise_qualimap_multiqc      misc                                     All default qualimap output for process :ref:`multiqc<Multiqc>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
\*                                       misc                                     All default qualimap output.             
======================================== ======================================== ========================================

**script**::

	qualimap bamqc -bam ${snippy_bam} --skip-duplicated -c -outformat "HTML" -outdir . -nt ${task.cpus}
	qualimapDir=${snippy_bam.baseName}_stats
	mv \$qualimapDir ${snippy_bam.baseName}

Nextstrain Metadata
-------------------

**script**::

	# The set command is to deal with PS1 errors
	set +eu
	# Enable conda activate support in this bash subshell
	CONDA_BASE=\$(conda info --base) ;
	source \$CONDA_BASE/etc/profile.d/conda.sh
	
	# Activate the nextstrain environment
	conda activate nextstrain-8.0.0
	
	# Format metadata
	${params.scriptdir}/format_metadata_Assembly.sh . ${sqlite} ${params.scriptdir}
	
	# Geocode
	divisions="country state"
	for div in \$divisions;
	do
	${params.scriptdir}/geocode_NextStrain.py \
	--in-tsv nextstrain/metadata_nextstrain.tsv \
	--loc-col BioSampleGeographicLocation \
	--out-tsv nextstrain/metadata_nextstrain_geocode_\${div}.tsv\
	--out-lat-lon nextstrain/lat_longs_\${div}.tsv \
	--div \${div};
	done
	
	cat \
	nextstrain/lat_longs_country.tsv \
	nextstrain/lat_longs_state.tsv > nextstrain/lat_longs_all.tsv
	
	
	# Deactivate the nextstrain environment
	conda deactivate

Nextstrain Treetime
-------------------

**script**::

	# The set command is to deal with PS1 errors
	set +eu
	# Enable conda activate support in this bash subshell
	CONDA_BASE=\$(conda info --base) ;
	source \$CONDA_BASE/etc/profile.d/conda.sh
	
	# Activate the nextstrain environment
	conda activate nextstrain-8.0.0
	
	mkdir -p nextstrain/treetime_clock/;
	treetime \
	--aln ${snippy_filter_aln} \
	--tree ${iqtree_treefile} \
	--dates ${metadata_nextstrain} \
	--clock-filter 3 \
	--keep-root \
	--gtr infer \
	--confidence \
	--keep-polytomies \
	--relax 1.0 0 \
	--max-iter 3 \
	--coalescent skyline \
	--covariation \
	--outdir nextstrain/treetime_clock \
	--date-column BioSampleCollectionDate \
	--verbose 6 2>&1 | tee nextstrain/treetime_clock/treetime_clock.log;
	
	# Deactivate env
	conda deactivate

Nextstrain Mugration
--------------------

**script**::

	# The set command is to deal with PS1 errors
	set +eu
	# Enable conda activate support in this bash subshell
	CONDA_BASE=\$(conda info --base) ;
	source \$CONDA_BASE/etc/profile.d/conda.sh
	
	# Activate the nextstrain environment
	conda activate nextstrain-8.0.0
	
	mkdir -p nextstrain/treetime_mugration_biovar/;
	mkdir -p nextstrain/treetime_mugration_country/;
	mkdir -p nextstrain/treetime_mugration_state/;
	
	treetime mugration \
	--tree ${timetree} \
	--attribute BioSampleBiovar \
	--states ${geocode_state} \
	--confidence \
	--outdir nextstrain/treetime_mugration_biovar/ \
	--verbose 6 2>&1 | tee nextstrain/treetime_mugration_biovar/treetime_mugration_biovar.log
	mv nextstrain/treetime_mugration_biovar/annotated_tree.nexus nextstrain/treetime_mugration_biovar/annotated_tree_biovar.nexus;
	mv nextstrain/treetime_mugration_biovar/confidence.csv nextstrain/treetime_mugration_biovar/confidence_biovar.csv  ;
	
	treetime mugration \
	--tree ${timetree} \
	--attribute country \
	--states ${geocode_state} \
	--confidence \
	--outdir nextstrain/treetime_mugration_country/ \
	--verbose 6 2>&1 | tee nextstrain/treetime_mugration_country/treetime_mugration_country.log
	mv nextstrain/treetime_mugration_country/annotated_tree.nexus nextstrain/treetime_mugration_country/annotated_tree_country.nexus;
	mv nextstrain/treetime_mugration_country/confidence.csv nextstrain/treetime_mugration_country/confidence_country.csv  ;
	
	treetime mugration \
	--tree ${timetree} \
	--attribute state \
	--states ${geocode_state} \
	--confidence \
	--outdir nextstrain/treetime_mugration_state/ \
	--verbose 6 2>&1 | tee nextstrain/treetime_mugration_state/treetime_mugration_state.log
	mv nextstrain/treetime_mugration_state/annotated_tree.nexus nextstrain/treetime_mugration_state/annotated_tree_state.nexus;
	mv nextstrain/treetime_mugration_state/confidence.csv nextstrain/treetime_mugration_state/confidence_state.csv  ;
	
	
	# Deactivate env
	conda deactivate

Nextstrain Json
---------------

**script**::

	# The set command is to deal with PS1 errors
	set +eu
	# Enable conda activate support in this bash subshell
	CONDA_BASE=\$(conda info --base) ;
	source \$CONDA_BASE/etc/profile.d/conda.sh
	
	# Activate the nextstrain environment
	conda activate nextstrain-8.0.0
	
	mkdir -p nextstrain/augur/;
	mkdir -p nextstrain/auspice/;
	
	augur refine \
	--alignment ${snippy_filter_aln} \
	--tree ${divergencetree} \
	--metadata ${metadata_nextstrain} \
	--output-tree nextstrain/augur/augur-refine.nwk \
	--output-node-data nextstrain/augur/mutation_lengths.json \
	--keep-root
	
	sed -i 's/branch_length/mutation_length/g' nextstrain/augur/mutation_lengths.json
	
	augur ancestral \
	--tree nextstrain/augur/augur-refine.nwk \
	--alignment ${snippy_core_vcf}  \
	--vcf-reference ${ref_chrom_fna} \
	--output-node-data nextstrain/augur/nt_muts.json \
	--output-vcf nextstrain/augur/augur-ancestral.vcf
	
	augur translate \
	--tree nextstrain/augur/augur-refine.nwk \
	--vcf-reference ${ref_chrom_fna} \
	--ancestral-sequences nextstrain/augur/augur-ancestral.vcf \
	--genes ${baseDir}/auspice/config/genes.txt \
	--reference-sequence ${ref_gff} \
	--output-node-data nextstrain/augur/aa_muts.json
	
	augur clades \
	--tree nextstrain/augur/augur-refine.nwk \
	--mutations nextstrain/augur/nt_muts.json \
	nextstrain/augur/aa_muts.json \
	--clades ${baseDir}/auspice/config/clades.csv \
	--output-node-data nextstrain/augur/clades.json
	
	${params.scriptdir}/treetime_dates_json.py \
	--time ${timetree} \
	--dates ${timetree_dates} \
	--json nextstrain/augur/branch_lengths.json
	
	${params.scriptdir}/treetime_mugration_json.py \
	--tree ${biovar_nexus} \
	--json nextstrain/augur/traits_biovar.json \
	--conf ${biovar_conf} \
	--trait biovar
	
	${params.scriptdir}/treetime_mugration_json.py \
	--tree ${country_nexus} \
	--json nextstrain/augur/traits_country.json \
	--conf  ${country_conf} \
	--trait country
	
	${params.scriptdir}/treetime_mugration_json.py \
	--tree ${state_nexus} \
	--json nextstrain/augur/traits_state.json \
	--conf ${state_conf} \
	--trait state
	
	augur export v2 \
	--tree nextstrain/augur/augur-refine.nwk \
	--metadata ${geocode_state} \
	--node-data nextstrain/augur/nt_muts.json \
	nextstrain/augur/aa_muts.json \
	nextstrain/augur/clades.json \
	nextstrain/augur/mutation_lengths.json \
	nextstrain/augur/branch_lengths.json \
	nextstrain/augur/traits_biovar.json \
	nextstrain/augur/traits_country.json \
	nextstrain/augur/traits_state.json \
	--output nextstrain/auspice/auspice.json \
	--lat-long ${lat_longs} \
	--auspice-config ${baseDir}/auspice/config/modernAssembly_auspice_config.json
	
	
	# Deactivate env
	conda deactivate

Multiqc
-------

Generate a MultiQC report from pipeline analyses.

======================================== ======================================== ========================================
Input:                                   Type                                     Description                              
======================================== ======================================== ========================================
ch_snippy_pairwise_qualimap_multiqc      misc                                     All default qualimap output from process :ref:`qualimap_snippy_pairwise<Qualimap_Snippy_Pairwise>`
======================================== ======================================== ========================================

======================================== ======================================== ========================================
Publish:                                 Type                                     Description                              
======================================== ======================================== ========================================
multiqc_report.html                      html                                     MultiQC report file.                     
*_data                                   misc                                     All default MultiQC data files.          
======================================== ======================================== ========================================

**script**::

	multiqc --config ${params.multiqc_config} .
