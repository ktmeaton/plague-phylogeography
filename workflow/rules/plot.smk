include: "functions.smk"
import os

rule plot_missing_data:
  input:
    filter_snp_log = expand(results_dir + "/snippy_multi/{{reads_origin}}/{{locus_name}}/filter{missing_data}/full/snippy-multi.snps.log",
		                         missing_data=config["snippy_multi_plot_missing_data"],
                           ),
  output:
    plot     = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/snippy-multi.snps.missing-data.html",
  log:
    notebook = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/processed_{locus_name}_plot_missing_data.py.ipynb",
  notebook:
	  os.path.join(notebooks_dir, "plot_missing_data.py.ipynb")

rule plot_snp_matrix:
  input:
    aln  = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/filter{missing_data}/{prune}/snippy-multi.snps.aln",
  output:
    dist = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/filter{missing_data}/{prune}/snippy-multi.snps.dist",
    html = results_dir + "/snippy_multi/{reads_origin}/{locus_name}/filter{missing_data}/{prune}/snippy-multi.snps.dist.heatmap.html",
  shell:
    """
    snp-dists {input.aln} > {output.dist}
    python3 {scripts_dir}/plot_distance_matrix.py -i {output.dist} -o {output.html}
    """
