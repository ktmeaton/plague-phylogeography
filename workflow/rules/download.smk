include: "functions.smk"
import os
#import itertools # Chaining list of lists of file accessions

# -----------------------------------------------------------------------------#
#                                Data Download                                 #
# -----------------------------------------------------------------------------#

rule download_sra:
  """
  Download SRA fastq files.
  """
  message: "Downloading and dumping {wildcards.sample} fastq for BioSample {wildcards.sample_dir}."
  output:
    fastq = results_dir + "/data/sra/{sample_dir}/{sample}_1.fastq.gz",
    log   = results_dir + "/data/sra/{sample_dir}/{sample}.log"
  resources:
    cpus = 1,
    time_min = 360,
  shell:
    "{scripts_dir}/download_sra.sh \
        {project_dir} \
        {results_dir}/data/sra/ \
        {wildcards.sample_dir} \
        {wildcards.sample} 1> {output.log}"

rule download_assembly:
    """
    Download files from the NCBI ftp server.
    """
    message: "Downloading and decompressing {wildcards.reads_origin} sample {wildcards.sample}.{wildcards.ext}"
    output:
        file = results_dir + "/data/{reads_origin}/{sample}/{sample}.{ext}"
    wildcard_constraints:
        ext = "(fna|gbff|gff)",
	      reads_origin = "(reference|assembly)",
    params:
		    ftp = lambda wildcards: globals()["identify_" + wildcards.reads_origin + "_ftp"]()
    resources:
        cpus = 1,
		shell:
        """
        for url in {params.ftp};
				do
				    if [[ $url =~ {wildcards.sample} ]]; then
						  sample_url=`echo $url | sed "s/fna/{wildcards.ext}/g"`;
							wget --quiet -O - $sample_url | gunzip -c > {output.file};
							if [[ {wildcards.reads_origin} == "reference" && ({wildcards.ext} == "fna" || {wildcards.ext} == "gff") ]]; then
							    python {scripts_dir}/rename_headers.py --file {output.file};
							fi
						fi;
				done
				"""

rule locus_bed:
  """
  Create a bed file of loci.
  """
  input:
    gbff = results_dir + "/data/{reads_origin}/{sample}/{sample}.gbff"
  output:
    bed = results_dir + "/data/{reads_origin}/{sample}/{sample}.bed"
  shell:
    """
    {scripts_dir}/locus_bed.sh {input.gbff} {output.bed}
    """
