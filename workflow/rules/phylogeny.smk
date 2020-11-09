import random # generate a random seed for iqtree

# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
  """
  Construct a maximum likelihood phylogeny.
  """
    input:
        full_aln = results_dir + "/snippy_multi/{reads_origin}/snippy-core_{locus_name}.full.aln",
        snp_aln = results_dir + "/snippy_multi/{reads_origin}/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
    output:
        tree = report(results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}.treefile",
                      caption=os.path.join(report_dir,"iqtree.rst"),
                      category="Phylogenetics",
                      subcategory="IQTREE"),
        log = report(os.path.join(logs_dir,"iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}.log"),
                     caption=os.path.join(report_dir, "logs.rst"),
                     category="Logs",
                     subcategory="Phylogeny"),
    params:
        #seed = random.randint(0, 99999999),
        seed = config["iqtree_seed"],
        prefix = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}",
    resources:
        load=100,
        time_min=600,
      	cpus=workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb=workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        os.path.join(logs_dir, "iqtree","{reads_origin}","iqtree-core_{locus_name}.filter" + "{missing_data}" + ".log")
    shell:
        "iqtree \
            -s {input.snp_aln} \
						{config[iqtree_model]} \
            --threads-max {resources.cpus} \
            -nt {resources.cpus} \
            -o {config[iqtree_outgroup]} \
            -seed {params.seed} \
            --runs {config[iqtree_runs]} \
            -fconst `snp-sites -C {input.full_aln}` \
            {config[iqtree_other]} \
            -pre {params.prefix} 1>{log}"
