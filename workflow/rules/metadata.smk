# Metadata filtering

rule metadata:
    """
		Create a metadata file of filename and strain name.
  	"""
		output:
		    tsv = results_dir + "/metadata/{reads_origin}/metadata.tsv",
		params:
		    samples = lambda wildcards: ",".join([
				      os.path.basename(os.path.dirname(path)) for path in identify_paths(outdir="metadata", reads_origin=wildcards.reads_origin)
						]),
				db = os.path.join(results_dir, "sqlite_db", config["sqlite_db"])
		shell:
		    """
		    python {scripts_dir}/metadata_file_strain.py \
				  --db {params.db} \
					--samples-csv {params.samples} \
					--output {output.tsv} ;
				"""
