import random # generate a random seed for iqtree

# -----------------------------------------------------------------------------#
#                         Maximum Likelihood Phylogeny                         #
# -----------------------------------------------------------------------------#
rule iqtree:
    """
    Construct a maximum likelihood phylogeny.
    """
    input:
        constant_sites = results_dir + "/snippy_multi/{reads_origin}/snippy-core_{locus_name}.full.constant_sites.txt",
        snp_aln        = results_dir + "/snippy_multi/{reads_origin}/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
    output:
        tree           = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}.treefile",
        iqtree         = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}.iqtree",
        log            = logs_dir    + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}.log",
    params:
        #seed = random.randint(0, 99999999),
        seed           = config["iqtree_seed"],
        prefix         = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}",
    resources:
        load           = 100,
        time_min       = 600,
      	cpus           = workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb         = workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        os.path.join(logs_dir, "iqtree","{reads_origin}","iqtree-core_{locus_name}.filter" + "{missing_data}" + ".log")
    shell:
        """
        iqtree \
            -s {input.snp_aln} \
						{config[iqtree_model]} \
            --threads-max {resources.cpus} \
            -nt {resources.cpus} \
            -o {config[iqtree_outgroup]} \
            -seed {params.seed} \
            --runs {config[iqtree_runs]} \
            -fconst `cat {input.constant_sites}` \
            {config[iqtree_other]} \
            -pre {params.prefix} 1>{log};
        """

rule iqtree_scf:
    """
    Estimate site concordance factors.
    """
    input:
        tree           = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}.treefile",
        snp_aln        = results_dir + "/snippy_multi/{reads_origin}/snippy-core_{locus_name}.snps.filter{missing_data}.aln",
        constant_sites = results_dir + "/snippy_multi/{reads_origin}/snippy-core_{locus_name}.full.constant_sites.txt",
    output:
        cf_branch      = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}_post.cf.branch",
        cf_stat        = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}_post.cf.stat",
        cf_tree        = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}_post.cf.tree",
    params:
        seed           = config["iqtree_seed"],
        prefix         = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}_post",
    resources:
        load           = 100,
        time_min       = 600,
      	cpus           = workflow.global_resources["cpus"] if ("cpus" in workflow.global_resources) else 1,
        mem_mb         = workflow.global_resources["mem_mb"] if ("mem_mb" in workflow.global_resources) else 4000,
    log:
        os.path.join(logs_dir, "iqtree","{reads_origin}","iqtree-core_{locus_name}.filter" + "{missing_data}_post.cf" + ".log")
    shell:
        """
        iqtree   \
            -t {input.tree}   \
            -s {input.snp_aln}   \
            --prefix {params.prefix}   \
            -fconst `cat {input.constant_sites}` \
            --scf 1000   \
            --seed {params.seed}   \
            --threads-max {resources.cpus} \
            -nt {resources.cpus};
        """

rule parse_tree:
    """
    Parse IQTREE trees and rename internal nodes.
    """
    input:
        tree     = results_dir + "/iqtree/{reads_origin}/iqtree-core_{locus_name}.filter{missing_data}_post.cf.tree",
        tsv      = results_dir + "/metadata/{reads_origin}/metadata.tsv",
    output:
        tsv      = results_dir + "/parse_tree/{reads_origin}/{locus_name}_filter{missing_data}/parse_tree.tsv",
        tree     = results_dir + "/parse_tree/{reads_origin}/{locus_name}_filter{missing_data}/parse_tree.nwk",
    log:
        notebook = results_dir + "/parse_tree/{reads_origin}/{locus_name}_filter{missing_data}/parse_tree_processed.py.ipynb",
    notebook:
        os.path.join(notebooks_dir, "parse_tree.py.ipynb")

rule clock_model:
    """
    Estimate a clock model from the parsed IQTREE phylogeny.
    """
    input:
        tree     = results_dir + "/parse_tree/{reads_origin}/{locus_name}_filter{missing_data}/parse_tree.nwk",
        tsv      = results_dir + "/parse_tree/{reads_origin}/{locus_name}_filter{missing_data}/parse_tree.tsv",
    output:
        tree     = results_dir + "/clock/{reads_origin}/{locus_name}_filter{missing_data}/clock_model_timetree.nwk",
        tsv      = results_dir + "/clock/{reads_origin}/{locus_name}_filter{missing_data}/clock_model.tsv",
    log:
        notebook = results_dir + "/clock/{reads_origin}/{locus_name}_filter{missing_data}/clock_model_processed.py.ipynb",
    notebook:
        os.path.join(notebooks_dir, "clock_model.py.ipynb")
