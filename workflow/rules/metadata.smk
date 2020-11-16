# Metadata filtering
include: "functions.smk"

rule metadata:
    """
		Create a metadata file of filename and strain name.
  	"""
		output:
		    txt = report(results_dir + "/metadata/{reads_origin}/metadata.txt",
                     caption=os.path.join(report_dir,"metadata.rst"),
                     category="Metadata",
                     subcategory="General"),
		params:
		    samples = lambda wildcards: ",".join(remove_duplicates([
				      os.path.basename(os.path.dirname(path)) for path in identify_paths(outdir="metadata", reads_origin=wildcards.reads_origin)
						])),
				db = os.path.join(results_dir, "sqlite_db", config["sqlite_db"])
		shell:
		    """
		    python {scripts_dir}/metadata.py \
				  --db {params.db} \
					--samples-csv {params.samples} \
					--output {output.txt} ;
				"""
