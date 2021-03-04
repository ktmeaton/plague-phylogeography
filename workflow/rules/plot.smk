include: "functions.smk"
import os

rule plot_missing_data:
  input:
    filter_snp_log = expand(logs_dir + "/snippy_multi/{{reads_origin}}/snippy-core_{{locus_name}}.snps.filter{missing_data}.log",
		                         missing_data=config["snippy_multi_plot_missing_data"],
                           ),
  output:
    plot = results_dir + "/snippy_multi/{reads_origin}/missing_data_{locus_name}.snps.html",
  log:
    notebook= os.path.join(logs_dir, "notebooks", "{reads_origin}",	"processed_{locus_name}_plot_missing_data.py.ipynb")
  notebook:
	  os.path.join(notebooks_dir, "plot_missing_data.py.ipynb")
